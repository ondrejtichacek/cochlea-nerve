function [ ANTStatistics, results, runopt ] = ANT( stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt, plotopt, args)
arguments
    stimulus
    topt
    midopt
    mechopt
    mnaopt
    antopt
    hhopt
    runopt
    opt
    memopt
    paropt
    plotopt

    args.init_flag = false
end
%ANT
disp_title( 'Auditory Nerve Transduction Simulation', '*' );

%% Init
disp_title( 'Synapse INIT' );

if args.init_flag
    
    [topt_ss_ihc_pot, stimulus_ss_ihc_pot] = devopts.stimulus( ...
        antopt.samplingFrequency, Time(100, 'ms'), ...
        'onset', Time(0, 'ms'), ...
        'offset', Time(0, 'ms'), ...
        'special', 'zero');

    [potential, displacement, TimeVoltage, Signal, TimeSignal, antopt, runopt] = ...
        iniANT( stimulus_ss_ihc_pot, topt_ss_ihc_pot, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt );

    topt.compute = topt.total;

    % use steady state
    potential = potential(end) * ones(size(stimulus.signal))';
    TimeVoltage = stimulus.time(:);
    Signal = stimulus.signal;
    TimeSignal = stimulus.time;

else

    [potential, displacement, TimeVoltage, Signal, TimeSignal, antopt, runopt] = ...
        iniANT( stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt );
end


xpos = antopt.positionsToExcite;

%%

%% Preprocessing
% * interpolate the potential (Voltage) Time and Signal
%   to respect the required sampling
% * scale voltage by provided factor

% disp('Using rand modification of IHC potential')
% potential = Vm_rand(:);

disp_title( 'Preprocessing' );

[potential, time, Signal] = ANTPreprocess( potential, TimeVoltage, Signal, TimeSignal, antopt, topt );

%% Synapse IC

use_IC = true;

if use_IC && ~args.init_flag
    extra_dirs = {runopt.read_only_cache};

    antopt_ic = copy(antopt);
    antopt_ic.numberOfRepetitions = 1;
    antopt_ic.slices = {antopt.slicesToExcite}; % THIS IS IMPORTANT
            
    [runopt, ~] = findResultExtended(extra_dirs, 'synapse_ic', ...
            [], [], [], mechopt, mnaopt, antopt_ic, [], runopt, opt);
    
    cache_ic_file = fullfile(runopt.path.synapse_ic, 'synapse_ic.mat');
    cache_ss_file = fullfile(runopt.path.synapse_ic, 'synapse_ss.mat');
    if ~exist(runopt.path.synapse_ic, 'dir')
        mkdir(runopt.path.synapse_ic)
    end
    
    if ~runopt.recalculate.synapse_ic ...
            && runopt.found.synapse_ic ...
            && exist(cache_ic_file, 'file') ...
            && exist(cache_ss_file, 'file')
        dbprintf('Loading initial conditions from a file:\n -- %s.\n', cache_ic_file);
        dbprintf('Loading initial conditions from a file:\n -- %s.\n', cache_ss_file);
        IC = load(cache_ic_file);
        SS = load(cache_ss_file);
    elseif runopt.no_create_only_load_synapse_ic
        
        error('We require initial conditions to exist.')
        
    else

        tf_init = Time(5000, 'ms');
        % tf_init = Time(1000, 'ms');

        tf_analysis_start = Time(4000, 'ms');
        % tf_analysis_start = Time(800, 'ms');
    
        [topt_ic, stimulus_ic] = devopts.stimulus( ...
                antopt.samplingFrequency, tf_init, ...
                'zeroDuration', Time(0, 'ms'), ...
                'fadeDuration', Time(0, 'ms'), ... 
                'onset', Time(0, 'ms'), ...
                'offset', Time(0, 'ms'), ...
                'special', 'zero');

        runopt_ic = copy(runopt);

        % runopt_ic.recalculate.ant_postprocess = true;

        runopt_ic.do.nerve = false;
        runopt_ic.plot.synapse = false;
        runopt_ic.plot.synapse_avg = false;
        runopt_ic.plot.nerve = false;

        runopt_ic.analysisStart = tf_analysis_start;
    
        if ~isempty(plotopt)
            plotopt_cpy = copy(plotopt);
        else
            plotopt_cpy = [];
        end

        [ANTStatistics_IC, res_IC] = ANT(stimulus_ic, topt_ic, ...
            copy(midopt), copy(mechopt), copy(mnaopt), antopt_ic, copy(hhopt), ...
            runopt_ic, copy(opt), copy(memopt), copy(paropt), plotopt_cpy, ...
            ...'xpos', xpos, ...
            'init_flag', true);
        f = load(ANTStatistics_IC.full.Properties.Source);
        
        tmp = [f.FreePool.TimeAverage];
        q0_mean = [tmp.MEAN];
        % q0 = round([tmp.MEAN]);
        
        tmp = [f.CaConc.TimeAverage];
        C0_mean = [tmp.MEAN];

        tmp = [f.ReStore.TimeAverage];
        w0_mean = [tmp.MEAN];
        % w0 = round([tmp.MEAN]);

        m0_mean = f.deterministic_stats.m.statistics.tspan.mean;
        
        IC = struct( ...
            'q0_mean', q0_mean, ...
            'C0_mean', C0_mean, ...
            'w0_mean', w0_mean, ...
            'm0_mean', m0_mean);

        t0_ic = find(res_IC.synapse.t > tf_analysis_start, 1);
        tf_ic = numel(res_IC.synapse.t);
        
        SS = struct( ...
            'q', res_IC.synapse.q(t0_ic:tf_ic,':'), ...
            'C', res_IC.synapse.C(t0_ic:tf_ic,':'), ...
            'w', res_IC.synapse.w(t0_ic:tf_ic,':'), ...
            'm', res_IC.synapse.m(t0_ic:tf_ic,':'));


        save(cache_ic_file, '-struct', 'IC', '-v7.3');
        save(cache_ss_file, '-struct', 'SS', '-v7.3');

        % Create SUCCESS indicator
        fid = fopen(fullfile(runopt.path.synapse_ic, 'SUCCESS'),'w');
        fclose(fid);
    end
    

else
    % IC for calculating IC
    IC = [];
    SS = [];
end

antopt.initialConditions = IC;
antopt.initialConditions_ss = SS;


%% Sometimes it is useful to finish just with IC calculation (serial) and do the rest in parallel

if runopt.calculate_only_ant_IC && ~args.init_flag
    ANTStatistics = [];
    results = [];
    return
end

%% Synapse
disp_title( 'Synapse Simulation' );

% get synapse result
[SynResFile, runopt, antopt] = handleSynapseResult( potential, time, stimulus, topt, midopt, mechopt, mnaopt, antopt, runopt, paropt, opt, memopt );

% save settings
% runopt = saveANTSettings( hhopt, runopt );


%% Nerve
if runopt.do.nerve
    disp_title( 'Nerve Simulation' );

    % get nerve result
    [NerveResFile, runopt, hhopt] = handleNerveResult( SynResFile, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, paropt, opt, memopt );
else
    NerveResFile = {};
end

%%

for i = 1:numel(SynResFile)
    SynResFile{i}.xgrid_ind = antopt.slicesToExcite;
    SynResFile{i}.xgrid = antopt.positionsToExcite;
    if runopt.do.nerve
        NerveResFile{i}.xgrid_ind = antopt.slicesToExcite;
        NerveResFile{i}.xgrid = antopt.positionsToExcite;
    end
end

%% postprocess result

% postprocess_windows = {'full', 'main', 'onset_20ms'};
% postprocess_windows = {'full', 'main'};

if isa(stimulus, 'PureTone')
    postprocess_windows = {...
        'full', 'main', ...
        'onset_10ms', 'offset_10ms', ...
        'onset_20ms', 'offset_20ms', ...
        'ss_zero_last_10ms', ...
        'ss_zero', 'ss_end_50ms'};
else
    postprocess_windows = {'full'};
end

if args.init_flag == true
    postprocess_windows = {'full'};
end


if runopt.do.ant_postprocess
    if numel(SynResFile{1}.n) > 1
        warning("Can't run analysis for multi-position sim.")
        ANTStatistics = [];
    else
        for i = 1:numel(postprocess_windows)
            w = postprocess_windows{i};
            StatFileName = fullfile(runopt.path.synapse, ...
                sprintf('ANTStatistics_%s.mat', w));
            
            doANTPostprocess = false;
            
            if runopt.recalculate.ant_postprocess || runopt.recalculate.synapse
                doANTPostprocess = true;
            end
        
            if ~exist(StatFileName,'file')
                doANTPostprocess = true;
            end
        
            if exist(StatFileName,'file')
                load(StatFileName, 'num_replicas_synapse', 'num_replicas_nerve');
                if ~exist('num_replicas_synapse', 'var') || ~exist('num_replicas_nerve', 'var')
                    warning('Variable num_replicas_synapse or num_replicas_nerve not found in ANTStatistics.mat, recalculating');
                    doANTPostprocess = true;
                else
                    if num_replicas_synapse < antopt.numberOfRepetitions
                        doANTPostprocess = true;
                    end
                    if num_replicas_nerve < antopt.numberOfRepetitions
                        doANTPostprocess = true;
                    end
                end
            end
            
            if doANTPostprocess
                ANTStatistics.(w) = ANTPostprocess( ...
                    SynResFile, NerveResFile, potential, time, Signal, ...
                    antopt, hhopt, runopt, stimulus, memopt, ...
                    'time_span', w);
            else
                ANTStatistics.(w) = matfile(StatFileName,'Writable', true);
            end
        end
    end
else
    ANTStatistics = [];
end


%% Plotting

if runopt.plot.ant_postprocess
    if isempty(ANTStatistics)
        warning('Can''t plot ANT Postprocess -- variable empty. Make sure runopt.do.ant_postprocess is set to true.')
    else
        for i = 1:numel(postprocess_windows)
            w = postprocess_windows{i};
            ANTPostprocessPlotting( ANTStatistics.(w), SynResFile, NerveResFile, potential, time, Signal, antopt, hhopt, runopt, plotopt, stimulus, 1 )
        end
    end
end

if runopt.plot.ant || runopt.plot.synapse || runopt.plot.synapse_avg || runopt.plot.nerve
    
    if ischar(antopt.plotSynapseSimulation)
        if strcmp(antopt.plotSynapseSimulation, 'all')
            L = 1 : antopt.numberOfRepetitions;
        end
    elseif isnumeric(antopt.plotSynapseSimulation)
        L = antopt.plotSynapseSimulation;
    else
        L = [];
    end
       
    SEL = [];
    for i = L
        if any(strcmp(SynResFile{i}.SR, runopt.SRsToPlot)) || any(strcmp('all', runopt.SRsToPlot))
            SEL(end+1) = i;
        end
    end
end

% plot ANT
if runopt.plot.ant 
    if numel(antopt.slicesToExcite) > 1
        ANTplotting( SynResFile(SEL), NerveResFile(SEL) )
    end
end

maxslice_relative = slicesToExcite({'max'}, mechopt.Numstacks, SynResFile{1}.V, [], mnaopt.xgrid);

% plot synapse results
if runopt.plot.synapse_avg
    disp_title( 'Synapse AVG Plotting' );
    fprintf('Plotting slice %d ...\n', antopt.slicesToExcite(maxslice_relative))

    if numel(SynResFile{1}.n) > 1
        warning("Can't plot for multi-position sim.")
    else
        for i = 1:numel(postprocess_windows)
            w = postprocess_windows{i};
            if isempty(NerveResFile)
                Nrv = [];
            else
                Nrv = NerveResFile(SEL);
            end
            SynapsePlotting_avg( potential, time, ANTStatistics.(w), SynResFile(SEL), Nrv, antopt, hhopt, opt, runopt, plotopt, maxslice_relative );
        end
    end
end
% plot synapse results
if runopt.plot.synapse
    disp_title( 'Synapse Plotting' );
    if numel(SEL) > 5
        SEL = SEL(1:5);
        warning('Attempting to plot too many simulations, the number reduced to 5')
    end
    
    if numel(SynResFile{1}.n) > 1
        warning("Can't plot for multi-position sim.")
    else
        fprintf('Plotting slice %d ...\n', antopt.slicesToExcite(maxslice_relative))
        for i = 1:numel(postprocess_windows)
            w = postprocess_windows{i};
            SynapsePlotting( potential, time, ANTStatistics.(w), SynResFile(SEL), NerveResFile(SEL), antopt, hhopt, opt, runopt, plotopt, maxslice_relative );
        end
    end
end

% plot nerve results
if runopt.plot.nerve
    disp_title( 'Nerve Plotting' );
    if numel(SEL) > 5
        SEL = SEL(1:5);
        warning('Attempting to plot too many simulations, the number reduced to 5')
    end
    if numel(NerveResFile{1}.n) > 1
        warning("Can't plot for multi-position sim.")
    else
        for i = 1:numel(SEL)
            nerveplotting(NerveResFile{SEL(i)}, stimulus, topt, mechopt, mnaopt, antopt, hhopt, runopt, opt, plotopt, memopt);
        end
    end
end

results.synapse = [SynResFile{:}];
results.nerve = [NerveResFile{:}];

end
