function [ ANTStatistics ] = ANTPostprocess( SynResFile, NerveResFile, ...
        Voltage, time, Signal, antopt, hhopt, runopt, stimulus, memopt, args)
%ANTPOSTPROCESS 
arguments
    SynResFile
    NerveResFile
    Voltage
    time
    Signal
    antopt
    hhopt
    runopt
    stimulus
    memopt
    args.time_span (1,:) char = 'full'
end

if runopt.analysisStart >= time(end)
    warning(['Can''t perform postprocessing (analysis) because reqeqsted analysis start\n', ...
    '\trunopt.analysisStart = %s\nis higher than computed\n\ttime span = [%s, %s]'], ...
        runopt.analysisStart, time(1), time(end));
    
    ANTStatistics = [];
    return
end

%% INIT

ANALYSIS = struct();

[Signal_Toolbox_licenceAvailable, ~] = license('checkout','Signal_Toolbox');

if ~Signal_Toolbox_licenceAvailable
    warning('Signal toolbox not available, skipping some parts of analysis!');
end

SR = getActiveSR(antopt);

nbins = 49;
if isa(stimulus, 'PureTone')
    nbins = min( floor( antopt.samplingFrequency/stimulus.frequency )+1 , nbins);
end

N = nbins + 1;

% analysisStartTime = runopt.analysisStart;
% analysisStartTimeIndex = min([length(time), find(analysisStartTime <= time,1)]);
% warning('Analysis domain striped to < %s , %s >.', time(analysisStartTimeIndex), time(end));

switch args.time_span
    case 'full'
        analysis_span = [time(1), time(end)];
        
    case 'main'
        analysis_span = stimulus.main_tspan;

    case 'onset_10ms'
        analysis_span = [stimulus.zeroDuration, stimulus.zeroDuration + Time(10, 'ms')];

    case 'offset_10ms'
        analysis_span = [stimulus.main_tspan(2), stimulus.main_tspan(2) + Time(10, 'ms')];

    case 'onset_20ms'
        analysis_span = [stimulus.zeroDuration, stimulus.zeroDuration + Time(20, 'ms')];

    case 'offset_20ms'
        analysis_span = [stimulus.main_tspan(2), stimulus.main_tspan(2) + Time(20, 'ms')];

    case 'ss_zero'
        analysis_span = [time(1), stimulus.zeroDuration];

    case 'ss_zero_last_10ms'
        analysis_span = [...
            max([Time(0), stimulus.zeroDuration - Time(10, 'ms')]), stimulus.zeroDuration];

    case 'ss_end_50ms'
        analysis_span = [stimulus.main_tspan(2) - Time(50, 'ms'), stimulus.main_tspan(2)];

    otherwise
        error('not implemented')
end

dbprintf('Analysing time span: %s\n', args.time_span)


%% Hilbert Transform
% * computes hilbert transform of Signal and returns
% * amplitude at Time
% * phase at Time
%

if Signal_Toolbox_licenceAvailable
    [~, SignalAmplitude, SignalPhase, InstantaneousFrequency] = HilbertTransform( Signal,  stimulus(1).samplingFrequency );
else
    warning('Can''t use hilbert transform without signal processing toolbox. Skipping this analysis.')
    [SignalAmplitude, SignalPhase] = deal([]);
end

%% Deterministic signal statistics

variables_to_analyse = struct( ...
        'V', struct( ...
            'description', 'membrane potential'), ...
        'C', struct( ...
            'description', 'Calcium concentration in the vicinity of the synaptic membrane'), ...
        'm', struct( ...
            'description', 'Fraction of open channels') ...
        );

ANALYSIS.deterministic_stats = time_statistics(SynResFile{1}, variables_to_analyse, runopt);

%%

%% Get average from replications -- SYNAPSE

% MinPeakDistance = min( 1/(2*MainSignalFrequency);
MinPeakDistance = eps;

PEAK_PROPERTIES_Synapse = { ...
    ...'NPeaks', 10000, ...
    ...'SortStr', 'ascend', ...
    ...'MinPeakHeight', 0.1, ...
    ...'MinPeakProminence', 0.1, ...
    ...'Threshold', 0.1, ...
    'MinPeakDistance', MinPeakDistance, ...
    ...'WidthReference', 'halfheight' ... 'halfprom'
    ...'MinPeakWidth', 0.1, ...
    ...'MaxPeakWidth', 0.1, ...
    };

F = griddedInterpolant(time.ms, SignalPhase);

% ttt = time;
ttt = Time(SynResFile{1,1}.ft.t, SynResFile{1,1}.ft.unit);

SignalPhase_synapse = F(ttt.ms);

analysis_span_ind(1) = find(analysis_span(1) <= ttt + Time(1e3*eps, 'ms'), 1);
analysis_span_ind(2) = find(analysis_span(2) <= ttt + Time(1e3*eps, 'ms'), 1);

dbprintf('Analysis domain striped to < %s , %s >.\n', ttt(analysis_span_ind(1)), ttt(analysis_span_ind(2)));

vars = {'c', 'q', 'w', 'k', 'C'};
vars = {'c_dyn', 'q_dyn', 'w_dyn', 'k', 'C'};
var_names = {'Cleft', 'FreePool', 'ReStore', 'NTRate', 'CaConc'};
do_find_peaks = logical([0, 0, 0, 0, 0]);



[AA, stat_time] = ReplicationAverage( ...
    SynResFile, ttt, PEAK_PROPERTIES_Synapse, ...
    runopt, antopt, hhopt, memopt, SignalPhase_synapse, N, ...
    analysis_span_ind, vars, var_names, ...
    'do_find_peaks', do_find_peaks); % no need to look for peaks

for i = 1:numel(vars)
    vname = var_names{i};
    ANALYSIS.(vname) = AA.(vname);
end

%% EPSC

PeakFields = PeakStatistics([], [], []);

stat = init_statistics(1, 1, SR, PeakFields, stat_time, 1);

it = 0;
for i = 1:numel(SynResFile)
    for j = 1:SynResFile{i}.n_dyn_rep
        it = it + 1;
        
        e = SynResFile{i}.vesicle_release_events_dyn{j};
        % e = SynResFile{i}.vesicle_release_events;
        % assert(numel(e) == 1)
        % e = e{1};
        e = e(:);
        % e = e(e > runopt.analysisStart.s);
        e = e(e >= analysis_span(1).s & e <= analysis_span(2).s);
        
        pks = ones(size(e));
        locs = e;
        
        [STAT, InterTimes] = PeakStatistics(stat_time, pks, locs);
    
        PeakFields = fields(STAT);
    
        stat.InterTimes{it} = InterTimes;
        stat.Peak{it} = pks;
        stat.Location{it} = locs;
        % stat.Width{it} = width;
        % stat.Prominence{it} = prom;
    
    
        for f = 1:length(PeakFields)
    
            [ stat.(PeakFields{f}).N, ...
              stat.(PeakFields{f}).MEAN, ...
              stat.(PeakFields{f}).M2, ...
              stat.(PeakFields{f}).VAR, ...
              stat.(PeakFields{f}).SD ] = OnlineVariance( STAT.(PeakFields{f}), ...
                                                       stat.(PeakFields{f}).N, ...
                                                       stat.(PeakFields{f}).MEAN, ...
                                                       stat.(PeakFields{f}).M2 );
            [ stat.(PeakFields{f}).MIN, ...
              stat.(PeakFields{f}).MAX ] = OnlineMinMax( STAT.(PeakFields{f}), ...
                                                       stat.(PeakFields{f}).MIN, ...
                                                       stat.(PeakFields{f}).MAX );
        end
    end
end

stat.unit = 's';
          
ANALYSIS.EPSC = stat;

%% Get average from replications -- NERVE

% MinPeakDistance = min( 1/(2*MainSignalFrequency);
% MinPeakDistance = eps;
MinPeakDistance = 0.5*1e-3; % = 0.5 ms

PEAK_PROPERTIES_Nerve = { ...
    ...'NPeaks', 10000, ...
    ...'SortStr', 'ascend', ...
    'MinPeakHeight', 80, ...
    ...'MinPeakProminence', 0.1, ...
    ...'Threshold', 0.1, ...
    'MinPeakDistance', MinPeakDistance, ...
    ...'WidthReference', 'halfheight' ... 'halfprom'
    ...'MinPeakWidth', 0.1, ...
    ...'MaxPeakWidth', 0.1, ...
    };

% Currenlty we don't calculate statistics on N_Na and Istim - needs to be
% reimplemented.
% vars = {'V', 'N_Na', 'Istim'};
% var_names = {'Volt', 'N_Na', 'Istim'};
vars = {'V'};
var_names = {'Volt'};
do_find_peaks = logical([1, 0, 0]);

if ~isempty(NerveResFile)
    
    F = griddedInterpolant(time.ms, SignalPhase);

    % ttt = time;
    ttt = Time(NerveResFile{1,1}.ft.t, NerveResFile{1,1}.ft.unit);
    
    SignalPhase_nerve = F(ttt.ms);

    analysis_span_ind(1) = find(analysis_span(1) <= ttt + Time(1e3*eps, 'ms'),1);
    analysis_span_ind(2) = find(analysis_span(2) <= ttt + Time(1e3*eps, 'ms'),1);
    
    dbprintf('Analysis domain striped to < %s , %s >.\n', ttt(analysis_span_ind(1)), ttt(analysis_span_ind(2)));

    [AA, stat_time] = ReplicationAverage( ...
        NerveResFile, ttt, PEAK_PROPERTIES_Nerve, ...
        runopt, antopt, hhopt, memopt, SignalPhase_nerve, N, ...
        analysis_span_ind, vars, var_names, ...
        'do_find_peaks', do_find_peaks); % no need to look for peaks


    for i = 1:numel(vars)
        vname = var_names{i};
        ANALYSIS.(vname) = AA.(vname);
    end

else
    for i = 1:numel(vars)
        vname = var_names{i};
        ANALYSIS.(vname) = [];
    end
end

%% Time-based statistics

% analysisStartTimeIndex = min([length(t),find(runopt.analysisStart.ms <= t,1)]);
% warning('Analysis domain striped to < %f , %f >.', t(analysisStartTimeIndex), t(end));
% 
% for i = 1:length(SR)
%     for j = 1:numel(stimulus)
%         [ TT_Cleft{i,j}, C_Cleft{i,j}, J_Cleft{i,j}, U_Cleft{i,j}, I_Cleft{i,j} ] = TimeIndexPhaseStat( Cleft.(SR{i}).MEAN, t, Voltage, analysisStartTimeIndex, stimulus(j) );    
%         if ~isempty(NerveResFile)
%             [ TT_Volt{i,j}, C_Volt{i,j}, J_Volt{i,j}, U_Volt{i,j}, I_Volt{i,j} ] = TimeIndexPhaseStat( Volt.(SR{i}).MEAN, t, Voltage, analysisStartTimeIndex, stimulus(j) );
%         end
%     end
% end


%% SAVE SELECTED VARIABLES TO FILE

StatFileName = fullfile(runopt.path.synapse, ...
   sprintf('ANTStatistics_%s.mat', args.time_span));

ANALYSIS.time = time;
ANALYSIS.PEAK_PROPERTIES_Synapse = PEAK_PROPERTIES_Synapse;
ANALYSIS.PEAK_PROPERTIES_Nerve = PEAK_PROPERTIES_Nerve;
ANALYSIS.Signal = Signal;
ANALYSIS.SignalAmplitude = SignalAmplitude;
ANALYSIS.SignalPhase = SignalPhase;
ANALYSIS.InstantaneousFrequency = InstantaneousFrequency;
ANALYSIS.analysis_span = analysis_span;
ANALYSIS.analysis_span_ind = analysis_span_ind;

save(StatFileName, '-struct', 'ANALYSIS', '-v7.3')

ANTStatistics = matfile(StatFileName,'Writable',true);

ANTStatistics.num_replicas_synapse = numel(SynResFile);
ANTStatistics.num_replicas_nerve = numel(NerveResFile);



end



