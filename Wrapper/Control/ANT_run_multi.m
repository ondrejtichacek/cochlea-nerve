function [analysis_variables, CACHE] = ANT_run_multi(varargin)
arguments (Repeating)
    varargin
end
arguments
end

args = struct();

[args, opt, memopt, paropt] = common_opts(args, varargin{:});

%% Plotting options

plotopt = plotOpt('MNA', ...
            ...
            'default_plt_action', false, ...
            ...
            ...'subplot', 0);
            'subplot', true);
        
% plotopt.hide_figures_while_plotting();

%% Default options

[~, stimulus, midopt, mechopt, mnaopt] = devopts.defaults(args, 'fun', plotopt);

%% Run options
runopt = runOpt( ...
    'CodeVersion', codeVersion(), ...
    'calculate_only_ant_IC', args.calculate_only_ant_IC, ...
    'Debug', false, ...
    'ReconstructResults', true, ...true, ...
    'draw_gui', false, ...
    'save_figures', false, ...
    'figureSaveDir', fullfile(opt.cochleadir, 'Results', 'ANT'), ...
    ...'SRsToPlot', {{'H1'}}, ... % must be as double cell {{}}
    'SRsToPlot', {{'all'}}, ... % must be as double cell {{}}
    'waitbarFunctionAvailable', false, ... % disable waitbar
    'verbose', 2 );

% verbose levels:   0 ... just matlab errors & warnings
%                   1 ... low verbosity
%                   2 ... normal verbosity
%                   3 ... high verbosity

runopt.no_create_only_load_oc_mna = args.no_create_only_load_oc_mna;
runopt.no_create_only_load = args.no_create_only_load;

% default action for do is TRUE
if isempty(args.skip_do)
    skip_do = { ...
    ...'mech', ...
    ...'oc_mna', ...
    ...'synapse', ...
    ...'nerve', ...
    ...'ant_postprocess' ...
    };
else
    skip_do = args.skip_do;
end

% default action for recalculate is FALSE
if isempty(args.do_recalculate)
    do_recalculate = { ...
        ...'mech', ...
        ...'mech_statistics', ...
        ...'oc_mna', ...
        ...'oc_mna_statistics', ...
        ...'synapse', ...
        ...'nerve', ...
        ...'replications', ...
        ...'ant_postprocess', ...
        };
else
    do_recalculate = args.do_recalculate;
end

% default action for plot is FALSE
if isempty(args.do_plot)
    do_plot = { ...
        ...'mech', ...
        ...'oc_mna', ...
        ...'ant', ...
        ...'synapse', ...
        ...'nerve', ...
        ...'ant_postprocess' ...
        };
else
    do_plot = args.do_plot;
end

% update in runopt
runopt.update_struct('do', skip_do, false);
runopt.update_struct('recalculate', do_recalculate, true);
runopt.update_struct('plot', do_plot, true);


%% Plotting options
% plotopt = plotOpt('LCRB', ...
%     ...
%     'JoinReplications', false, ...
%     ...
%     'subsamplingFactor', 1, ... 10
%     'subsamplingFactorFilt', 1, ... 50
%     'subsamplingFactorStores', 1, ... 50
%     'movingAverageWindowLength', 0.01, ... [s]
%     'tspan', [0, tf] ...
%     );

%% Plotting options
plotopt = plotOpt('LCRB', ...
    ...
    'default_plt_action', false, ...
    ...
    'JoinReplications', false, ...
    ...
    'subsamplingFactor', 20, ... 10
    'subsamplingFactorFilt', 100, ... 50
    'subsamplingFactorStores', 100, ... 50
    'movingAverageWindowLength', Time(10, 'ms'), ...
    ...'tspan', [Time(0), tf] ...
    'tspan', [Time(0, 'ms'), Time(100, 'ms')] ...
    );

% default action for plotopt.do is FALSE
do_draw = { ...
    ... 'Channels', ...
    ... 'Vars1', ...
    ... 'Vars2', ...
    ... 'Stores', ...
    ... 'Cleft', ...
    ... 'CleftNerve', ...
    ...
    ... 'nerve_maximal_crossection', ...
    ... 'nerve_voltage_mesh', ...
    ...
    ...'MNAsettings' ...
    ...
    ...'FourierTransform', ...
    ...'MaximalCrossection', ...
    ...
    ...'CurrentMesh', ...
    ...'VoltageMesh', ...
    ...
    ...'CurrentSteadyState', ...
    ...'VoltageSteadyState', ...
    ...
    ...'VoltageMesh_HC', ...
    ...'VoltageSteadyState_HC', ...
    ...
    ...'BMdispl', ...
    ...'BMdisplProfile', ...
    ...'TMdispl', ...
    ...'TMdisplProfile', ...
    ...
    ...'stimulus', ...
    };

% update in plotopt
plotopt.update_struct('do', do_draw, true);

plotopt.hide_figures_while_plotting();

%% ANT

ANT_fun = 'ANT';
% ANT_fun = 'ANT_single';
% ANT_fun = 'ANT_clamp';

ANTSamplingFrequency = 0.5 * args.GlobalSamplingFrequency;
HHSamplingFrequency = max(Frequency(100, 'kHz'), args.GlobalSamplingFrequency);

slices = args.slices;
if isempty(args.slices)
    slices = {{'max_std_new'}};
end

    function antopt = set_antopt(CONF)

        antopt = antOpt( ...
            'NumDiv', 1, ...
            'save_method', 'matlab_matfile', ...
            'script', 'ANT', ...
            ...
            ...'ant', 'eguia', ...
            ...'ant', 'meddis', ...
            ...'ant', 'sumner', ...
            ...'ant', 'sumnerStochastic', ...
            ...'ant', 'dev', ...
            ...'ant', 'dev_stochastic', ...
            'ant', args.ant_ver, ...
            'fiber', 'default', ...
            ...
            'plotSynapseSimulation', false, ... % false, index vector (e.g. 1:3) or 'all'
            'plotSynapseVarsAreIdentical', true, ... % true if plotting replication of just one synapse parameter set
            ...
            'samplingFrequency', ANTSamplingFrequency, ... sampling frequency of the Synapse result & Nerve input and result
            ...
            ...'force_quick_onset', true', ...
            ...
            ... % ! needs to be two "cell" {} brackets as struct function "removes one" to extend dimensions !!!
            ...
            ...'slices', {{round(linspace(1,mnaopt.Numstacks,100))}}, ...
            ...'slices', {{'all'}}, ... % all slices
            ...'slices', {{'max'}}, ... % slice with maximal Curr.IHC value
            ...'slices', {{'max+',3}}, ... % max-slice + right Neighbourhood
            ...'slices', {{'max-',3}}, ... % max-slice + left Neighbourhood
            ...'slices', {{'max+-',1}}, ... % max-slice + Neighbourhood
            ...'slices', {{25:75}}, ... % multiple slices direct selection
            ...'slices', {{256}}, ... % direct slice selection
            ...'slices', {{'max_std'}}, ...
            'slices', slices, ...
            ...
            'VoltageAmplitudeFactor', 1, ...
            'VoltageAmplitudeAddFactor', 0, ...
            'IgnoreNTRateThreshold', false, ...
            ...
            'numberOfRepetitions', args.num_replicas ...
            );

       if regexp(CONF.fiber_properties{1}.SR, 'ring_variable_(\d+)')
            antopt.transductionopt = transductionOpt_v4_1(antopt.ant, CONF.fiber_properties{1});
            antopt.fiber = antopt.transductionopt.SR;

        else
            antopt.transductionopt = transductionOpt_v4(antopt.ant, CONF.fiber_properties{1}.SR);
            antopt.fiber = antopt.transductionopt.SR;
        end
    end

hhopt = hhOpt( ...
    'NumDiv', 1, ...
    'save_method', 'matlab_matfile', ...
    'method', args.hh_method, ...
    ...'method', 'hhnint', ...  % deterministic HH
    ...'method', 'hhfloor', ...  % deterministic HH
    ...'method', 'fxnint', ...  % stochastic HH - Fox algorithm
    ...'method', 'cw', ...      % stochastic HH - Chow&White algorithm
    ...
    'script', 'ANT', ...
    ...
    'plotNerveSimulation', false, ... % false, index vector (e.g. 1:3) or 'all'
    ...
    'samplingFrequency', HHSamplingFrequency, ... sampling frequency of the Synapse result & Nerve input and result
    ...
    'scaleFactor', 2*1e-5 ... % sumnerstoch
    ...'scaleFactor', v2struct(1e-4, 1e-4, 1e-4, fieldNames) ... % sumnerstoch
    );

%%

do = struct( ...
    'analysis', 1, ...
    'gathering', 1, ...
    'rate_intensity',  0, ...
    'rate_intensity_2', 0, ...
    'rate_intensity_3', 0, ...
    'inter_times_distribution',  0, ...
    'inter_times_heatmap',  0, ...
    'phase_distribution',  0, ... (!) computationally intensive
    'hist_over_signal',  0, ...
    'hist_over_signal_masking', 0, ...
    'hist_over_signal_masking_2', 0, ...
    'hist_over_signal_masking_4', 0, ...
    'hist_over_signal_inc_dec', 0, ...
    'fano_factor',  0, ...
    'inp_signal_phase',  0 ...
    );

if args.cached_analysis == true
    do.analysis = false;
    % do.gathering = false;
end
if isempty(args.cache)
    CACHE = struct();
else
    CACHE = args.cache;
end

% SRS = args.SRS;

% for i = 1:numel(args.fiber_properties)
%     SRS{i} = args.fiber_properties{i}.SR;
% end

SRS = [];



%% ANALYSIS

UPDATE_OPTS = args.update_opts;

% update opts (eg when comparing two setups)
if ~isempty(UPDATE_OPTS)
    for i = 1:2:numel(UPDATE_OPTS)
        opt_name = UPDATE_OPTS{i};
        opt_val = UPDATE_OPTS{i+1};
        switch opt_name
            case 'do'
                do = opt_val;
            case 'mnaopt'
                mnaopt = opt_val;
            otherwise
                error('opt name %s not recognised', opt_name)
        end
    end
end

% override_analysis_start = Time(0);
override_analysis_start = [];

keys = fields(plotopt.do);
for i = 1:numel(keys)
    key = keys{i};
    plotopt.do.(key) = false;
end

PLOT_FUN = @ANT_plot;
PLOT_FUN = [];

[Configurations, Conf_nums, Conf_info] = ParameterProduct( ...
        args.product_args{:});

if isfield(Configurations, 'fiber_properties')
    for i = 1:numel(Configurations)
        fiber_dislay_name{i} = Configurations(i).fiber_properties{1}.DisplayName;
        Configurations(i).fiber_display_name = fiber_dislay_name{i};
    end
end

ConfTable = struct2table(Configurations);

%%

AMPLITUDE = [];
if ~isempty(args.Amplitude)
    AMPLITUDE = args.Amplitude;
elseif isfield(Conf_info.range, 'amplitude')
    AMPLITUDE = Conf_info.range.amplitude;
end

FREQUENCY = [];
if ~isempty(args.Frequency)
    FREQUENCY = args.Frequency;
elseif isfield(Conf_info.range, 'frequency')
    FREQUENCY = Conf_info.range.frequency;
end

% DRY_RUN = true;
DRY_RUN = false;


%% Plot stimulus 

if false
    figure
    hold on 
    for i = 1:numel(Configurations)

        Conf = Configurations(i);
        [~, sig] = stimulus(Conf);

        % envelope = sig.envelope_fcn('ms');
        % plot(sig.audiotime.ms, envelope(sig.audiotime.ms) * 10^(sig.amplitude/20))
        
        % plot(sig.audio * 10^(sig.amplitude/20))

        envelope = sig.eval_fcn('ms');
        plot(sig.audiotime.ms, envelope(sig.audiotime.ms) * 10^(sig.amplitude/20))
    end
end

%%

if do.analysis

    SIM_OPTS = {midopt, mechopt, mnaopt, @set_antopt, hhopt};
    
    stats = batchGathering(...
        @(varargin) getANTPostprocessStatistics(ANT_fun, varargin{:}), ...
                1, Configurations, ...
                stimulus, SIM_OPTS, ...
                runopt, opt, memopt, paropt, ...
                'skip_configuration_uniqueness_test', true, ...
                'override_analysis_start', override_analysis_start, ...
                'fig_save_dir_name_fun', @(conf) fullfile(opt.cochleadir, 'Results', 'ANT'), ...
                'parfeval_iteration_delay', args.parfeval_iteration_delay, ...
                'PLOT_FUN', PLOT_FUN, ...
                'plotopt', plotopt);
       
    variables = cell(Conf_nums);
    for i = 1:prod(Conf_nums)
        % sub = Conf_info.sub(i,:);
        variables{i} = stats{i};
    end

    CACHE.variables = variables;
else
    variables = CACHE.variables;
end

analysis_variables = [];

% return

%% Gather Results

% SR = getActiveSR(antopt);
SR = SRS;

assert(numel(args.analysis_plot_variables) == 1)
args.analysis_plot_variables = args.analysis_plot_variables{1};

switch args.analysis_plot_variables
    case 'none'
        return
    case 'NTRate'

        SIGNAL = {'NTRate'};
        SIGNAL_TO_PLOT = {'NTRate'};
        VARS_s = {'TimeAverage'};
        VARS_q = {'phase'};
        VARS_d = {};
        VARS_t = {'time', 'SignalPhase', 'analysis_span', 'analysis_span_ind'};
        STATS = {'N', 'MEAN', 'SD', 'MAX', 'MIN', 'VAR'};
        postprocess_window = 'main';

        plot_scale = 1;
        unit = 'Hz';

        % do.rate_intensity =  1;
        do.rate_intensity_2 = 1;
        % do.inter_times_distribution =  1;
        do.inter_times_heatmap =  1;
        do.phase_distribution =  1;
        % do.hist_over_signal =  1;
        % do.inp_signal_phase =  1;

    case 'AP'

        SIGNAL = {'Volt'};
        SIGNAL_TO_PLOT = {'Volt'};
        VARS_s = {'TimeAverage'};
        VARS_q = {};
        VARS_d = {'InterTimes', 'Location'};
        % VARS_t = {'time', 'SignalPhase', 'analysis_span', 'analysis_span_ind'};
        % STATS = {'N', 'MEAN', 'SD', 'MAX', 'MIN', 'VAR'};
        VARS_t = {};
        STATS = {'MEAN'};

        % postprocess_window = 'full';
        postprocess_window = args.postprocess_window;

        plot_scale = 1;
        unit = 'Hz';

        % do.rate_intensity =  1;
        % do.rate_intensity_2 = 1;
        % do.inter_times_distribution =  1;
        % do.inter_times_heatmap =  1;
        % do.phase_distribution =  1;

        if strcmp(postprocess_window, 'full')
            do.hist_over_signal =  1;
        end
%         if strcmp(postprocess_window, 'main')
            % do.fano_factor =  1;
%         end
        % do.inp_signal_phase =  1;


    case 'EPSC'

        SIGNAL = {'EPSC'};
        SIGNAL_TO_PLOT = {'EPSC'};
        VARS_s = {'GlobalPeakRate'};
        VARS_q = {};
        VARS_d = {'InterTimes', 'Location'};
        VARS_t = {'time', 'SignalPhase', 'analysis_span', 'analysis_span_ind'};
        STATS = {'N', 'MEAN', 'SD', 'MAX', 'MIN', 'VAR'};

        % postprocess_window = 'main';
        postprocess_window = args.postprocess_window;

        plot_scale = 1;
        unit = 'Hz';
        
        % do.rate_intensity =  1;
        % do.rate_intensity_2 = 1;
        % do.inter_times_distribution =  1;
        % do.inter_times_heatmap =  1;
        % do.phase_distribution =  1;
        
        if strcmp(postprocess_window, 'full')
            do.hist_over_signal =  1;
        end
%         if strcmp(postprocess_window, 'main')
            do.fano_factor =  1;
%         end

        % do.inp_signal_phase =  1;
    
    case 'EPSC_full'

        SIGNAL = {'EPSC'};
        SIGNAL_TO_PLOT = {'EPSC'};
        VARS_s = {'GlobalPeakRate'};
        VARS_q = {};
        VARS_d = {'InterTimes', 'Location'};
        VARS_t = {'time', 'SignalPhase', 'analysis_span', 'analysis_span_ind'};
        STATS = {'N', 'MEAN', 'SD', 'MAX', 'MIN', 'VAR'};

        postprocess_window = 'full';

        plot_scale = 1;
        unit = 'Hz';

        % do.rate_intensity =  1;
        do.rate_intensity_2 = 1;
        % do.inter_times_distribution =  1;
        % do.inter_times_heatmap =  1;
        % do.phase_distribution =  1;
        do.hist_over_signal =  1;
        % do.inp_signal_phase =  1;

    case 'EPSC_onset_20ms'

        SIGNAL = {'EPSC'};
        SIGNAL_TO_PLOT = {'EPSC'};
        VARS_s = {'GlobalPeakRate'};
        VARS_q = {};
        VARS_d = {'InterTimes', 'Location'};
        VARS_t = {'time', 'SignalPhase', 'analysis_span', 'analysis_span_ind'};
        STATS = {'N', 'MEAN', 'SD', 'MAX', 'MIN', 'VAR'};

        postprocess_window = 'onset_20ms';

        plot_scale = 1;
        unit = 'Hz';

        % do.rate_intensity =  1;
        do.rate_intensity_2 = 1;
        % do.inter_times_distribution =  1;
        % do.inter_times_heatmap =  1;
        % do.phase_distribution =  1;
        % do.hist_over_signal =  1;
        % do.inp_signal_phase =  1;

    case 'CaConc'
        
        SIGNAL = {'CaConc'};
        SIGNAL_TO_PLOT = {'CaConc'};
        VARS_s = {'TimeAverage'};
        VARS_q = {'phase'};
        VARS_d = {};
        VARS_t = {'time', 'SignalPhase', 'analysis_span', 'analysis_span_ind'};
        STATS = {'N', 'MEAN', 'SD', 'MAX', 'MIN', 'VAR'};
        postprocess_window = 'full';

        plot_scale = 1e6;
        unit = 'uM';

        do.hist_over_signal = false;

        % do.rate_intensity =  1;
        % do.rate_intensity_2 = 1;
        % do.inter_times_distribution =  1;
        % do.inter_times_heatmap =  1;
        % do.phase_distribution =  1;
        % do.hist_over_signal =  1;
        % do.inp_signal_phase =  1;

    otherwise
        error('Unknown analysis_plot_variables: %s', args.analysis_plot_variables)
end


if ~isempty(args.do_override)
    fn = fieldnames(args.do_override);
    for i = 1:numel(fn)
        do.(fn{i}) = args.do_override.(fn{i});
    end
end


% STATS = {};
% STATS = {'MEAN', 'SD'};
% SIGNAL = {'CaConc', 'Cleft'};
% SIGNAL = {'Cleft'};
% SIGNAL = {'Cleft', 'Volt'};
% SIGNAL = {'Cleft', 'NTRate'};
% SIGNAL = {'Cleft', 'Volt', 'NTRate'};
% SIGNAL = {'Cleft'};

% SIGNAL_TO_PLOT = {'Cleft'};
% SIGNAL_TO_PLOT = {'Cleft', 'Volt'};


% deterministic_signals = {'C'};

VARS_TO_LOAD = [SIGNAL, VARS_t{:}];

if do.gathering
    
    % prealoc
    t.gathering_prealoc = tic;
    for r = 1:numel(SIGNAL)
        for p = 1:numel(VARS_s)
            for q = 1:numel(STATS)
                RES.(SIGNAL{r}).(VARS_s{p}).(STATS{q}) = zeros(Conf_nums);
            end
        end
        for p = 1:numel(VARS_d)
            RES.(SIGNAL{r}).(VARS_d{p}) = cell(Conf_nums);
        end
        for p = 1:numel(VARS_t)
            RES.(VARS_t{p}) = cell(Conf_nums);
        end
    end
    fprintf('Gathering prealoc completed in %s\n',disp_toc(toc(t.gathering_prealoc)));
    
    deterministic_stats = cell(Conf_nums);
    
    t.gathering = tic;
    for i = 1:prod(Conf_nums)
        if isempty(variables{i})
            warning('Statistics for simulation %d   not found', i);
        else
            S = load(variables{i}.(postprocess_window).Properties.Source, VARS_TO_LOAD{:} );

            for r = 1:numel(SIGNAL)
                for p = 1:numel(VARS_s)
                    for q = 1:numel(STATS)
                        tmp = [S.(SIGNAL{r}).(VARS_s{p})];
                        tmp = [tmp.(STATS{q})];
                        % tmp = mean([tmp.(STATS{q})]);
                        if any(strcmp(SIGNAL{r}, {'CaConc'}))
                            tmp = mean(tmp);
                        elseif any(strcmp(SIGNAL{r}, {'NTRate'}))
                            tmp = sum(tmp);
                        end
                        RES.(SIGNAL{r}).(VARS_s{p}).(STATS{q})(i) = tmp; 

                    end
                end
                for p = 1:numel(VARS_q)
                    tmp = [S.(SIGNAL{r}).(VARS_q{p})];
                    tmp = [tmp.mean];
                    if any(strcmp(SIGNAL{r}, {'CaConc'}))
                        tmp = mean(tmp, 2);
                    elseif any(strcmp(SIGNAL{r}, {'NTRate'}))
                        tmp = mean(tmp, 2);
                    end
                    RES.(SIGNAL{r}).(VARS_q{p}){i} = tmp;
                end
                for p = 1:numel(VARS_d)
                    RES.(SIGNAL{r}).(VARS_d{p}){i} = cat(1,S.(SIGNAL{r}).(VARS_d{p}){:});
                    RES_2.(SIGNAL{r}).(VARS_d{p}){i} = {S.(SIGNAL{r}).(VARS_d{p}){:}};
                end
            end

            DET = load(variables{i}.(postprocess_window).Properties.Source, 'deterministic_stats');
            deterministic_stats{i} = DET.deterministic_stats;

            done_frac = i / prod(Conf_nums);
            time_expected = (1-done_frac)*toc(t.gathering)/(done_frac);
            fprintf('%.2f%% done, %s expected\n', 100*done_frac,disp_toc(time_expected));
        end
    end
    fprintf('Gathering completed in %s\n',disp_toc(toc(t.gathering)));

    for i = 1:prod(Conf_nums)
        if ~isempty(deterministic_stats{i})
            stats_templ = deterministic_stats{i};
            if isa(stats_templ, 'matlab.io.MatFile')
                stats_templ = load(stats_templ.Properties.Source);
            end
            stats_templ = nanstruct(stats_templ);
            break % from inner most loop
        end
    end
    for i = 1:prod(Conf_nums)
        if isempty(deterministic_stats{i})
            deterministic_stats{i} = stats_templ;
        end
    end

    %%
    % for i = 1:numel(AMPLITUDE)
    %     for j = 1:numel(FREQUENCY)
    %         for s = 1:numel(SRS)
    % 
    %             if isempty(variables{i,j,s})
    %                 warning('Statistics for simulation %d-%d-%d   not found', i, j, s);
    %             else
    %                 S = load(variables{i,j,s}.(postprocess_window).Properties.Source, VARS_TO_LOAD{:} );                
    % 
    %                 for r = 1:numel(SIGNAL)
    % %                 for s = 1:numel(SR)
    %                     for p = 1:numel(VARS_s)
    %                         for q = 1:numel(STATS)
    %                             tmp = [S.(SIGNAL{r}).(VARS_s{p})];
    %                             tmp = [tmp.(STATS{q})];
    %                             % tmp = mean([tmp.(STATS{q})]);
    %                             if any(strcmp(SIGNAL{r}, {'CaConc'}))
    %                                 tmp = mean(tmp);
    %                             elseif any(strcmp(SIGNAL{r}, {'NTRate'}))
    %                                 tmp = sum(tmp);
    %                             end
    %                             RES.(SIGNAL{r}).(SRS{s}).(VARS_s{p}).(STATS{q})(i,j) = tmp; 
    % 
    %                         end
    %                     end
    %                     for p = 1:numel(VARS_q)
    %                         tmp = [S.(SIGNAL{r}).(VARS_q{p})];
    %                         tmp = [tmp.mean];
    %                         if any(strcmp(SIGNAL{r}, {'CaConc'}))
    %                             tmp = mean(tmp, 2);
    %                         elseif any(strcmp(SIGNAL{r}, {'NTRate'}))
    %                             tmp = mean(tmp, 2);
    %                         end
    %                         RES.(SIGNAL{r}).(SRS{s}).(VARS_q{p}){i,j} = tmp;
    %                     end
    %                     for p = 1:numel(VARS_d)
    %                         RES.(SIGNAL{r}).(SRS{s}).(VARS_d{p}){i,j} = cat(1,S.(SIGNAL{r}).(VARS_d{p}){:});
    %                         RES_2.(SIGNAL{r}).(SRS{s}).(VARS_d{p}){i,j} = {S.(SIGNAL{r}).(VARS_d{p}){:}};
    %                     end
    % %                 end
    %                 end
    % 
    %                 for p = 1:numel(VARS_t)
    %                     RES.(VARS_t{p}){i,j} = S.(VARS_t{p});
    %                 end
    % 
    %                 DET = load(variables{i,j,s}.(postprocess_window).Properties.Source, 'deterministic_stats');
    %                 deterministic_stats{i,j,s} = DET.deterministic_stats;
    %             end
    %         end
    % 
    %         done_frac = (numel(FREQUENCY)*(i-1)+j)/numel(AMPLITUDE)/numel(FREQUENCY);
    %         time_expected = (1-done_frac)*toc(t.gathering)/(done_frac);
    %         fprintf('%.2f%% done, %s expected\n', 100*done_frac,disp_toc(time_expected));
    %     end
    % end
    % fprintf('Gathering completed in %s\n',disp_toc(toc(t.gathering)));
    
    % for i = 1:numel(AMPLITUDE)
    %     for j = 1:numel(FREQUENCY)
    %         for s = 1:numel(SRS)
    %             if ~isempty(deterministic_stats{i,j,s})
    %                 stats_templ = deterministic_stats{i,j,s};
    %                 if isa(stats_templ, 'matlab.io.MatFile')
    %                     stats_templ = load(stats_templ.Properties.Source);
    %                 end
    %                 stats_templ = nanstruct(stats_templ);
    %                 break % from inner most loop
    %             end
    %         end
    %     end
    % end
    % for i = 1:numel(AMPLITUDE)
    %     for j = 1:numel(FREQUENCY)
    %         for s = 1:numel(SRS)
    %             if isempty(deterministic_stats{i,j,s})
    %                 deterministic_stats{i,j,s} = stats_templ;
    %             end
    %         end
    %     end
    % end
end

%%

analysis_variables = RES;

if args.noplot
    dbprintf('Not plotting, finishing.')
    return
end
if ~canplot()
    dbprintf('Plotting not available, finishing.')
    return
end

%% Plots
if false
if strcmp(SIGNAL{r}, 'CaConc')
    for s = 1:numel(SRS)
        for i = 1:numel(AMPLITUDE)
            for j = 1:numel(FREQUENCY)
                t = RES.time{i,j};
                if isempty(variables{i,j,s})
                    warning('not found')
                else
                    S = load(variables{i,j,s}.(postprocess_window).Properties.Source, 'CaConc');
                    S = S.CaConc;
                    figure
                    hold on
                    for k = 1:numel(S)
                        plot(t.ms, S(k).MEAN * 1e6, ...
                            'Color', [0.5,0.5,0.5])
                    end
                    m = mean([S.MEAN], 2);
                    plot(t.ms, m * 1e6, ...
                        'Color', 'r')
                    xlabel('Time (ms)')
                    ylabel('Ca2+ (uM)')
                    title(sprintf('%s %d Hz %g dB SPL', SRS{s}, FREQUENCY(j), AMPLITUDE(i)))
                end
            end
        end
    end
end
end
%%

if do.rate_intensity_3
    
    postprocess_windows = { ...
        'main', ...
        'onset_10ms', ...
        'offset_10ms', ...
        ...'onset_20ms', 
        ...'offset_20ms', ...
        ...'ss_zero', 
        'ss_zero_last_10ms', ...
        'ss_end_50ms', ...
        ...'full', ...
        };

    SS = {};
    for ipw = 1:numel(postprocess_windows)
        pw = postprocess_windows{ipw};
        for i = 1:prod(Conf_nums)
            SS{i}.(pw) = load(variables{i}.(pw).Properties.Source, 'Volt', 'analysis_span'); %, VARS_TO_LOAD{:} );
        end
    end

    FIB = unique({Configurations.fiber_display_name});
    
    for r = 1:numel(SIGNAL_TO_PLOT)
    for i_fib = 1:numel(FIB)

        sel = strcmp(ConfTable.fiber_display_name, FIB{i_fib});

        CFT = ConfTable(sel,:);
        CFT = table2struct(CFT)';

        figure()
        hold on

        [~, cmap] = linecolors(numel(postprocess_windows));

        for ipw = 1:numel(postprocess_windows)
            pw = postprocess_windows{ipw};

            for i = 1:numel(CFT)
                Conf = CFT(i);
                [~, stim] = stimulus(Conf);
                A(i) = stim.amplitude;
    
                dt = diff(SS{i}.(pw).analysis_span);

                N_mean(i) = SS{i}.(pw).Volt.NumberOfPeaks.MEAN / dt.s;
                N_sd(i) = SS{i}.(pw).Volt.NumberOfPeaks.SD / dt.s;

                % N(i) = numel(cat(1, SS{i}.(pw).Volt.Location{:}));                
                % n_dyn_rep = 100;
                % N(i) = N(i) / args.num_replicas / n_dyn_rep;
            end

            plot_error = false;

            if ~plot_error
                plot(A, N_mean, ...
                    'Color', cmap(ipw,:), ...
                    'DisplayName', pw)
            else
                shadedErrorBar( ...
                    A,...
                    N_mean, N_sd / sqrt(6400), ...
                    'lineProps', ...
                        {'Color', cmap(ipw,:), ...
                        'DisplayName', pw, ...
                        'LineStyle', '-', ...
                        'LineWidth', 0.75}, ...
                    'transparent', true ...
                    )
            end
        end

        RES_tmp = RES;
        RES_tmp.(SIGNAL_TO_PLOT{r}).Location = RES_tmp.(SIGNAL_TO_PLOT{r}).Location(sel);

        data = analysis.ant.hist_over_signal_new(RES_tmp, SIGNAL_TO_PLOT{r}, ...
            CFT, Conf_nums, Conf_info,...
            runopt, opt, args.num_replicas, stimulus, ...
            'plotflag', false, ...
            'BinWidth', Time(5, 'ms'));

        M = cellfun(@(x) max(x.counts), data, 'UniformOutput', true);
        plot(A, M, ...
            'DisplayName', 'onset max', ...
            'Color', cmap(2,:))

        M = cellfun(@(x) mean(x.counts(end-3:end-1)), data, 'UniformOutput', true);
        plot(A, M, ...
            'DisplayName', 'steady state', ...
            'Color', cmap(4,:))

        xlabel('Level (dB SPL)')
        ylabel('Firing rate (spikes/s)')

        legend('location', 'best')
    end
    end
    
end


if do.rate_intensity
    if numel(AMPLITUDE) > 1
        analysis.ant.rate_intensity(RES, FREQUENCY, AMPLITUDE, SR, VARS_s, VARS_d, VARS_t, STATS, SIGNAL, runopt, opt)
    else
        warning('Can''t plot IO plot when numel(AMPLITUDE) = %d', numel(AMPLITUDE))
    end
end
if do.rate_intensity_2
    if numel(AMPLITUDE) > 1
        analysis.ant.rate_intensity_2(RES, FREQUENCY, AMPLITUDE, SR, VARS_s, SIGNAL, runopt, opt, args.amplitude_unit, ...
            "plot_scale", plot_scale);
    else
        warning('Can''t plot IO plot when numel(AMPLITUDE) = %d', numel(AMPLITUDE))
    end

end

if do.rate_intensity_2
    if numel(AMPLITUDE) > 1
        for s = 1:numel(SRS)
            Vmax = Voltage(cellfun(@(S) S.V.statistics.tspan.max, deterministic_stats(:,:,s)), deterministic_stats{1}.V.unit)';
            Vmin = Voltage(cellfun(@(S) S.V.statistics.tspan.min, deterministic_stats(:,:,s)), deterministic_stats{1}.V.unit)';
            V0 = Voltage(cellfun(@(S) S.V.statistics.init.val0, deterministic_stats(:,:,s)), deterministic_stats{1}.V.unit)';
        
            VV = Voltage(zeros(size(V0)));
            for i = 1:numel(VV)
                if abs(V0(i) - Vmax(i)) > abs(V0(i) - Vmin(i))
                    VV(i) = Vmax(i);
                else
                    VV(i) = Vmin(i);
                end
            end
            hfig = figure();
            yyaxis left
            analysis.ant.rate_intensity_2(RES, FREQUENCY, VV.mV, SR(s), VARS_s, SIGNAL, runopt, opt, args.amplitude_unit, ...
                "hfig", hfig, ...
                "plot_scale", plot_scale)
        
            yyaxis right
            V = linspace(-80, 10, 100);
            plot(V, boltzmann(V,1,-25.8, 5.5), 'r')
            plot(V, boltzmann(V,1,-31.4, 7.7), 'b')
        end
        
        hfig = figure;
        for i = 1:numel(FREQUENCY)
            KK = {i};
            analysis.ant.rate_intensity_simple(RES, FREQUENCY, VV.mV, SR, VARS_s, SIGNAL, ...
                runopt, opt, args.amplitude_unit, KK, ...
                "hfig", hfig, ...
                "plot_scale", plot_scale);
        end
        legend()
            xlim([min(VV(:).mV), max(VV(:).mV)])
            xlabel('Voltage step (mV)')
            ylabel('Release rate (Hz)')
    
        [Configurations, Conf_nums] = ParameterProduct( ...
            'Amplitude', AMPLITUDE, ...
            'Frequency', FREQUENCY, ...
            'fiber', SRS{1} );
    
        [topt, sig] = stimulus(Configurations(1));
    
        
        tt = sig.audiotime;
        fun = sig.envelope_fcn('ms');
        
        voltage_step = false;
        if voltage_step
            figure
            hold on
            for i = 1:numel(VV)
                plot(tt.ms, V0(i).mV + AMPLITUDE(i)* fun(tt.ms), 'k')
            end
            xlabel('Time (ms)')
            ylabel('Voltage (mV)')
        end
    else
        warning('Can''t plot IO plot when numel(AMPLITUDE) = %d', numel(AMPLITUDE))
    end
end

%%

% -------------------------------------------------------------------------
if false
signal = 'C'; % Ca concentration

figure
hold on
analysis.ant.var_intensity_profile( FREQUENCY, AMPLITUDE, deterministic_stats, signal )
title('intensity profile')
subtitle('peak Ca^{2+} concentration')
xlabel('Level [dB SPL]')
ylabel('Ca conc. max [\muM]')
end
% -------------------------------------------------------------------------
if false
if numel(AMPLITUDE) > 1
    signal = 'm'; % channels fraction
    
    figure
    hold on
    analysis.ant.var_intensity_profile( FREQUENCY, AMPLITUDE, deterministic_stats, signal )
    title('intensity profile')
    subtitle('fraction of open ?? channels')
    xlabel('Level [dB SPL]')
    ylabel('')
end
end
% -------------------------------------------------------------------------
if false
if numel(AMPLITUDE) > 1
    signal = 'V'; % voltage
    
    figure
    hold on
    analysis.ant.var_intensity_profile( FREQUENCY, AMPLITUDE, deterministic_stats, signal, 'mV' )
    title('intensity profile')
    subtitle('IHC basolateral voltage')
    xlabel('Level [dB SPL]')
    ylabel('Voltage. max [mV]')
end
end
% -------------------------------------------------------------------------

% xlim([0,20])

%%
if false
    signal_x = 'V'; % voltage
    signal_y = 'C'; %

    figure
    hold on
    analysis.ant.var_var_profile( FREQUENCY, AMPLITUDE, deterministic_stats, signal_x, signal_y, 'mV', '' )
    title('IHC voltage profile')
    subtitle('peak Ca^{2+} concentration')
    xlabel('Voltage. max [mV]')
    ylabel('Ca conc. max [\muM]')
end

%% INTER TIMES

if do.inter_times_distribution
    
%     SIGNAL_TO_PLOT = {'Cleft'};
%     SIGNAL_TO_PLOT = {'Cleft', 'Volt'};
%     AMPLITUDE_TO_PLOT = cellfun(@(x) find(x==AMPLITUDE), {90});
    AMPLITUDE_TO_PLOT = 1:numel(AMPLITUDE);
%     FREQUENCY_TO_PLOT = cellfun(@(x) find(x==FREQUENCY), {600 2000 10000});
    FREQUENCY_TO_PLOT = 1:numel(FREQUENCY);
%     SR_TO_PLOT = {'H1'};%{'H1', 'M1', 'L1'};
    SR_TO_PLOT = SR;
    
    if is_in('InterTimes', VARS_d)
        analysis.ant.inter_times_distribution(RES, FREQUENCY, AMPLITUDE, SIGNAL, ...
            FREQUENCY_TO_PLOT, AMPLITUDE_TO_PLOT, SR_TO_PLOT, SIGNAL_TO_PLOT, ...
            runopt, opt);
    end
    
end

if do.inter_times_heatmap
    
%     SIGNAL_TO_PLOT = {'Cleft'};
%     SIGNAL_TO_PLOT = {'Cleft', 'Volt'};
%     AMPLITUDE_TO_PLOT = cellfun(@(x) find(x==AMPLITUDE), {90});
    AMPLITUDE_TO_PLOT = 1:numel(AMPLITUDE);
%     FREQUENCY_TO_PLOT = cellfun(@(x) find(x==FREQUENCY), {600 2000 10000});
    FREQUENCY_TO_PLOT = 1:numel(FREQUENCY);
%     SR_TO_PLOT = {'H1'};%{'H1', 'M1', 'L1'};
    SR_TO_PLOT = SR;
    
    if is_in('InterTimes', VARS_d)
        analysis.ant.inter_times_heatmap(RES, FREQUENCY, AMPLITUDE, SIGNAL, ...
            FREQUENCY_TO_PLOT, AMPLITUDE_TO_PLOT, SR_TO_PLOT, SIGNAL_TO_PLOT, ...
            runopt, opt);
    end
    
end

%%

if do.phase_distribution

%     SIGNAL_TO_PLOT = {'Cleft'};
%     AMPLITUDE_TO_PLOT = cellfun(@(x) find(x==AMPLITUDE), {30,60,90});
    AMPLITUDE_TO_PLOT = 1:numel(AMPLITUDE);
%     SR_TO_PLOT = {'H1', 'M1', 'L1'};
%     SR_TO_PLOT = {'H1'};
    SR_TO_PLOT = SR;
    
    if is_in('Location', VARS_d)
        nbins = 49;
        analysis.ant.phase_distribution(RES, FREQUENCY, AMPLITUDE, ...
            AMPLITUDE_TO_PLOT, SR_TO_PLOT, SIGNAL_TO_PLOT, ...
            ANTSamplingFrequency, nbins, ...
            runopt, opt)
    end
    if is_in('phase', VARS_q)
        if isempty(args.special_signal)
            nbins = 49;
            analysis.ant.phase_average(RES, FREQUENCY, AMPLITUDE, ...
                AMPLITUDE_TO_PLOT, SR_TO_PLOT, SIGNAL_TO_PLOT, ...
                ANTSamplingFrequency, nbins, ...
                runopt, opt)
        end
    end
end

%%
if do.fano_factor
    if is_in('Location', VARS_d)
        for r = 1:numel(SIGNAL_TO_PLOT)
            
            hfig = figure();
            hold on
            set(gca, 'xscale', 'log')
            set(gca, 'yscale', 'log')
    
            [x, y] = Peterson_2014.fig4a();

            [x1, y1, x2, y2] = Bruce_2018.fig8(plotflag=false);

            plot(x, y, 'ko-', ...
                'DisplayName', 'experiment - Peterson (2014)')

            plot(x1, y1, 'Color', [0.75, 0.75, 0.75], ...
                'DisplayName', 'simulation - Zilany (2014)')

            plot(x2, y2, 'Color', [0.5, 0.5, 0.5], ...
                'DisplayName', 'simulation - Bruce (2018)')

            FIB = unique({Configurations.fiber_display_name});
            
            % cmap = lcmp(numel(FIB));
            [~, cmap] = linecolors(1);
        
            markers = {'o', 'x', 's', 'd', 'p'};
        
            for i_fib = 1:numel(FIB)
        
                sel = strcmp(ConfTable.fiber_display_name, FIB{i_fib});
        
                CFT = ConfTable(sel,:);
                CFT = table2struct(CFT)';
        
                RES_2_tmp = RES_2;
                Location = RES_2_tmp.(SIGNAL_TO_PLOT{r}).Location(sel);

                RES_tmp = RES;
                RES_tmp.(SIGNAL_TO_PLOT{r}).Location = RES_tmp.(SIGNAL_TO_PLOT{r}).Location(sel);
            
                for i = 1:numel(Location)
                    [T, F] = analysis.ant.fano_factor(Location{i});
                end

                figure(hfig)
                plot(T, F, ...
                    'Color', cmap(i_fib,:), ...
                    'LineWidth', 1, ...
                     ...'Marker', markers{i_fib}, ...
                    'DisplayName', sprintf('sim., fiber %s', FIB{i_fib}))
                
                % p = polyfit(log(T), log(F), 5);
                % x1 = logspace(0, 3, 100);
                % y1 = polyval(p, log(x1));
                % plot(x1, exp(y1));

                for i = 1:numel(Location)
                    [T, F] = analysis.ant.fano_factor(Location{i}, 'shuffle', true);
                end

                figure(hfig)
                plot(T, F, ...
                    'Color', cmap(i_fib,:), ...
                    'LineWidth', 1, ...
                    'LineStyle', ':', ...
                     ...'Marker', markers{i_fib}, ...
                    'DisplayName', sprintf('shuffled., fiber %s', FIB{i_fib}))

                if true
                    hfig2 = figure();
                    tl = tiledlayout('flow');
                    title(tl, FIB{i_fib});
                    
                    analysis.ant.hist_over_signal_new(RES_tmp, SIGNAL_TO_PLOT{r}, ...
                        CFT, Conf_nums, Conf_info,...
                        runopt, opt, args.num_replicas, stimulus, ...
                        'hfig', hfig2, ...
                        'plotflag', true, ...
                        'BinWidth', Time(1, 'ms'));
                end

            end
            
            figure(hfig)
            set(gca, 'Yscale', 'log')
            set(gca, 'Xscale', 'log')
            xlabel('Counting time T (ms)')
            ylabel('Fano factor')
            legend('Location', 'best')
            ylim([0.2,10])
        end
    end
end

%%
if do.hist_over_signal
    if is_in('Location', VARS_d)
        for r = 1:numel(SIGNAL_TO_PLOT)

                bw = Time(1, 'ms');

                if isempty(AMPLITUDE)
                    conf_DisplayName = [];
                else
                    conf_DisplayName = cellfun(@(txt) sprintf('%d dB SPL', txt), num2cell(AMPLITUDE), 'UniformOutput', false);
                end
    
                hist_data{r} = analysis.ant.hist_over_signal_new(RES, SIGNAL_TO_PLOT{r}, ...
                        Configurations, Conf_nums, Conf_info,...
                        runopt, opt, args.num_replicas, stimulus, ...
                        'conf_DisplayName', conf_DisplayName, ...
                        'BinWidth', bw);
   
                
                % for i = 1:numel(Configurations)
                %     Conf = Configurations(i);
                % end
        end
    end
end

%%
if do.hist_over_signal_masking
    if is_in('Location', VARS_d)
        for r = 1:numel(SIGNAL_TO_PLOT)

            hfig = figure;
            hold on
            set(gca, 'xscale', 'log')

            xlabel('\Delta T (ms)')
            ylabel('Probe response magnitude (% control)')

            [x, y] = Harris_1979.fig3(plotflag=false);
            plot(x, y, 'bo--', 'DisplayName', 'Harris and Dallos (1979), unit 34-20');

            [x, y, z] = Harris_1979.fig14a(plotflag=false);
            ind = find(z == 100, 1);
            plot(x, y(:,ind), 'bx--', 'DisplayName', 'Harris and Dallos (1979), unit 36-24');

            [x, y, z] = Harris_1979.fig14b(plotflag=false);
            ind = find(z == 100, 1);
            plot(x, y(:,ind), 'bs--', 'DisplayName', 'Harris and Dallos (1979), unit 39-1');

            [x, y, z] = Harris_1979.fig14c(plotflag=false);
            ind = find(z == 200, 1);
            plot(x, y(:,ind), 'bd--', 'DisplayName', 'Harris and Dallos (1979), unit 37-16');

            [x, y, z] = Harris_1979.fig9a(plotflag=false);
            ind = find(z == 30, 1);
            y = y(:,ind);
            plot(x(~isnan(y)), y(~isnan(y)), 'bp--', 'DisplayName', 'Harris and Dallos (1979), unit 40-7');

            [x, y, z] = Harris_1979.fig15a(plotflag=false);
            ind = find(z == 100, 1);
            plot(x, y(:,ind), 'rx-', 'DisplayName', 'Harris and Dallos (1979), median response');

            [x, y, z] = Harris_1979.fig10(plotflag=false);
            ind = find(z == 30, 1);
            % plot(x, y(:,ind), 'ro-', 'DisplayName', 'Harris and Dallos (1979), median response');

            ylim([0,105])
            xlim([0.9,330])

            legend('Location', 'SouthEast');


            FIB = unique({Configurations.fiber_display_name});
            
            cmap = lcmp(numel(FIB));

            markers = {'o', 'x', 's', 'd', 'p'};

            for i_fib = 1:numel(FIB)

                sel = strcmp(ConfTable.fiber_display_name, FIB{i_fib});

                CFT = ConfTable(sel,:);
                CFT = table2struct(CFT)';

                RES_tmp = RES;
                RES_tmp.(SIGNAL_TO_PLOT{r}).Location = RES_tmp.(SIGNAL_TO_PLOT{r}).Location(sel);

                M = zeros(1, numel(CFT));
    
                for i = 1:numel(CFT)
                    Conf = CFT(i);
                    Location = Time(RES_tmp.(SIGNAL_TO_PLOT{r}).Location{i}, 's');
                    [~, stim] = stimulus(Conf);
                    if Conf.delay < inf
                        tsp = stim.compounds(2).main_tspan;
                        % tsp = tsp - stim.compounds(2).onsetDuration; % see onset/offset in the paper
                    else
                        tsp = stim.main_tspan;
                    end

                    tsp = tsp + Time(2, 'ms'); % SYSTEM DELAY (?)
                    tsp(1) = tsp(1) - Time(1, 'ms'); % ONSET


                    M(i) = sum(Location.ms >= tsp(1).ms & Location.ms <= tsp(2).ms);
                end
    
                figure(hfig)

                plot(Conf_info.range.delay + 1, M/M(end)*100, ...
                    'Color', 'k', ...
                    'LineWidth', 1, ...
                    'Marker', markers{i_fib}, ...
                    'DisplayName', sprintf('simulation, fiber %s', FIB{i_fib}))
                
                if true
                    hfig2 = figure();
                    tl = tiledlayout('flow');
                    title(tl, FIB{i_fib});

                    delay = num2cell(Conf_info.range.delay + 1);
                    conf_DisplayName = cellfun(@(txt) sprintf('\\Delta = %d ms', txt), delay, 'UniformOutput', false);

                     analysis.ant.hist_over_signal_new(RES_tmp, SIGNAL_TO_PLOT{r}, ...
                        CFT, Conf_nums, Conf_info,...
                        runopt, opt, args.num_replicas, stimulus, ...
                        'conf_DisplayName', conf_DisplayName, ...
                        'hfig', hfig2, ...
                        'plotflag', true, ...
                        'BinWidth', Time(1, 'ms'));
                end

            end

            if false
                bws = Time(linspace(0.25,0.75,5), 'ms');
                for ii = 1:numel(bws)
    
                    bw = bws(ii);
        
                    hist_data{r} = analysis.ant.hist_over_signal_new(RES, SIGNAL_TO_PLOT{r}, ...
                            Configurations, Conf_nums, Conf_info,...
                            runopt, opt, args.num_replicas, stimulus, ...
                            'plotflag', false, ...
                            'BinWidth', bw);
                    
                    for i = 1:numel(Configurations)
                        Conf = Configurations(i);
                        [~, stim] = stimulus(Conf);
                        if Conf.delay < inf
                            tsp = stim.compounds(2).main_tspan;
                        else
                            tsp = stim.main_tspan;
                        end
                        i1 = find(hist_data{r}{i}.edges >= tsp(1).ms, 1);
                        i2 = find(hist_data{r}{i}.edges >= tsp(2).ms, 1);
        
                        tmp = max(hist_data{r}{i}.counts(i1:i2));
                        if isempty(tmp)
                            tmp = NaN;
                        end
                        M(i) = tmp;
                    end
                    
                    MMM{ii} = M/M(end)*100;
                end
    
                MMM = cat(1, MMM{:});
    
                figure(hfig)
                errorbar(Conf_info.range.delay + 1, mean(MMM,1), std(MMM,1), 'ko');
                errorbar(Conf_info.range.delay + 1, mean(MMMs,1), std(MMMs,1), 'kx');
            end
        end
    end
end

%%

if do.hist_over_signal_masking_2
    if is_in('Location', VARS_d)
        for r = 1:numel(SIGNAL_TO_PLOT)

            FIB = unique({Configurations.fiber_display_name});

            for i_fib = 1:numel(FIB)

                Harris_1979.fig10(plotflag=true, plot_fit=false)
                legend
                hold on
                hfig = gcf();
        
                THR = unique(ConfTable.masker_ampl_re_thr);
                
                cmap = lcmp(numel(THR));
        
                FIB = unique({Configurations.fiber_display_name});
        
                markers = {'o', 'x', 's', 'd', 'p'};

                for i_thr = 1:numel(THR)
    
                    sel = strcmp(ConfTable.fiber_display_name, FIB{i_fib}) & ...
                        ConfTable.masker_ampl_re_thr == THR(i_thr,:);
    
                    CFT = ConfTable(sel,:);
                    CFT = table2struct(CFT)';
                
                    RES_tmp = RES;
                    RES_tmp.(SIGNAL_TO_PLOT{r}).Location = RES_tmp.(SIGNAL_TO_PLOT{r}).Location(sel);
    
                    M = zeros(1, numel(CFT));
    
                    for i = 1:numel(CFT)
                        Conf = CFT(i);
                        Location = Time(RES_tmp.(SIGNAL_TO_PLOT{r}).Location{i}, 's');
                        [~, stim] = stimulus(Conf);
                        if Conf.delay < inf
                            tsp = stim.compounds(2).main_tspan;
                            % tsp = tsp - stim.compounds(2).onsetDuration; % see onset/offset in the paper
                        else
                            tsp = stim.main_tspan;
                        end

                        % tsp = tsp + Time(2, 'ms'); % SYSTEM DELAY (?)

                        M(i) = sum(Location.ms >= tsp(1).ms & Location.ms <= tsp(2).ms);
                    end
        
                    figure(hfig)
                    plot(Conf_info.range.delay + 1, M/M(end)*100, ...
                        'Color', cmap(i_thr,:), ...
                        'LineWidth', 1, ...
                        'Marker', markers{i_fib}, ...
                        'DisplayName', sprintf('simulation, fiber %s', FIB{i_fib}))
                    
                    if true

                        hfig2 = figure();
                        tl = tiledlayout('flow');
                        title(tl, FIB{i_fib});

                        analysis.ant.hist_over_signal_new(RES_tmp, SIGNAL_TO_PLOT{r}, ...
                            CFT, Conf_nums, Conf_info,...
                            runopt, opt, args.num_replicas, stimulus, ...
                            'hfig', hfig2, ...
                            'plotflag', true, ...
                            'BinWidth', Time(1, 'ms'));
                    end
    
                end
            end
        end
    end
end

%%

if do.hist_over_signal_masking_4
    if is_in('Location', VARS_d)
        for r = 1:numel(SIGNAL_TO_PLOT)
            
            Harris_1979.fig14a(plotflag=true)
            legend
            hold on
            hfig = gcf();

            DUR = unique(ConfTable.masker_duration);
            
            cmap = lcmp(numel(DUR));

            FIB = unique({Configurations.fiber_display_name});
            
            % cmap = lcmp(numel(FIB));

            markers = {'o', 'x', 's', 'd', 'p'};

            for i_fib = 1:numel(FIB)

                for i_dur = 1:numel(DUR)
    
                    sel = strcmp(ConfTable.fiber_display_name, FIB{i_fib}) & ...
                        ConfTable.masker_duration == DUR(i_dur,:);
    
                    CFT = ConfTable(sel,:);
                    CFT = table2struct(CFT)';
                
                    RES_tmp = RES;
                    RES_tmp.(SIGNAL_TO_PLOT{r}).Location = RES_tmp.(SIGNAL_TO_PLOT{r}).Location(sel);
    
                    M = zeros(1, numel(CFT));
    
                    for i = 1:numel(CFT)
                        Conf = CFT(i);
                        Location = Time(RES_tmp.(SIGNAL_TO_PLOT{r}).Location{i}, 's');
                        [~, stim] = stimulus(Conf);
                        if Conf.delay < inf
                            tsp = stim.compounds(2).main_tspan;
                            % tsp = tsp - stim.compounds(2).onsetDuration; % see onset/offset in the paper
                        else
                            tsp = stim.main_tspan;
                        end

                        % tsp = tsp + Time(2, 'ms'); % SYSTEM DELAY (?)

                        M(i) = sum(Location.ms >= tsp(1).ms & Location.ms <= tsp(2).ms);
                    end
        
                    figure(hfig)
                    plot(Conf_info.range.delay + 1, M/M(end)*100, ...
                        'Color', cmap(i_dur,:), ...
                        'LineWidth', 1, ...
                        'Marker', markers{i_fib}, ...
                        'DisplayName', sprintf('simulation, fiber %s', FIB{i_fib}))
                    
                    if true
                        hfig2 = figure();
                        tl = tiledlayout('flow');
                        title(tl, FIB{i_fib});
                        
                        analysis.ant.hist_over_signal_new(RES_tmp, SIGNAL_TO_PLOT{r}, ...
                            CFT, Conf_nums, Conf_info,...
                            runopt, opt, args.num_replicas, stimulus, ...
                            'hfig', hfig2, ...
                            'plotflag', true, ...
                            'BinWidth', Time(1, 'ms'));
                    end
    
                end
            end
        end
    end
end

%%

if do.hist_over_signal_inc_dec
    if is_in('Location', VARS_d)
        for r = 1:numel(SIGNAL_TO_PLOT)
        SD = linspace(1.5,2.5,20);
        for i_sd = 1:numel(SD)
            Zilany_2009.fig9b(plot=true)
            hfig = gcf();
            hold on


            FIB = unique({Configurations.fiber_display_name});
            
            cmap = lcmp(numel(FIB));

            markers = {'o', 'x', 's', 'd', 'p'};

            for i_fib = 1:numel(FIB)

                sel = strcmp(ConfTable.fiber_display_name, FIB{i_fib});

                CFT = ConfTable(sel,:);
                CFT = table2struct(CFT)';

                RES_tmp = RES;
                RES_tmp.(SIGNAL_TO_PLOT{r}).Location = RES_tmp.(SIGNAL_TO_PLOT{r}).Location(sel);

                M1 = zeros(1, numel(CFT));
                M2 = zeros(1, numel(CFT));
    
                for i = 1:numel(CFT)
                    Conf = CFT(i);
                    Location = Time(RES_tmp.(SIGNAL_TO_PLOT{r}).Location{i}, 's');
                    [~, stim] = stimulus(Conf);
                    if Conf.delay < inf
                        tsp = stim.compounds(2).main_tspan;
                        % tsp = tsp - stim.compounds(2).onsetDuration; % see onset/offset in the paper
                    else
                        tsp = stim.main_tspan;
                    end

                    % tsp = tsp + Time(2, 'ms'); % SYSTEM DELAY (?)
                    tsp = tsp + Time(SD(i_sd), 'ms'); % SYSTEM DELAY (?)
                    
                    tsp(1) = tsp(1) - Time(1, 'ms'); % ONSET
                    
                    tsp_1 = copy(tsp);
                    tsp_2 = copy(tsp);

                    tsp_1(2) = tsp_1(1) + Time(0.64, 'ms');
                    tsp_2(2) = tsp_2(1) + Time(10.2, 'ms');

                    M1(i) = sum(Location.ms >= tsp_1(1).ms & Location.ms <= tsp_1(2).ms);
                    T = tsp_1(2) - tsp_1(1);
                    M1(i) = M1(i) / 6400 / T.s;

                    M2(i) = sum(Location.ms >= tsp_2(1).ms & Location.ms <= tsp_2(2).ms);
                    T = tsp_2(2) - tsp_2(1);
                    M2(i) = M2(i) / 6400 / T.s;
                end
    
                figure(hfig)

                plot(Conf_info.range.delay, M1 - M1(end), ...
                    'Marker', 'o', ...
                    'Color', 'k', ...
                    'LineWidth', 2, ...
                    'DisplayName', sprintf('simulation, fiber %s', FIB{i_fib}))

                plot(Conf_info.range.delay, M2 - M2(end), ...
                    'Marker', 'x', ...
                    'Color', 'k', ...
                    'LineWidth', 2, ...
                    'DisplayName', sprintf('simulation, fiber %s', FIB{i_fib}))
                
                if false
                    hfig2 = figure();
                    tl = tiledlayout('flow');
                    title(tl, FIB{i_fib});

                    delay = num2cell(Conf_info.range.delay + 1);
                    conf_DisplayName = cellfun(@(txt) sprintf('\\Delta = %d ms', txt), delay, 'UniformOutput', false);

                     analysis.ant.hist_over_signal_new(RES_tmp, SIGNAL_TO_PLOT{r}, ...
                        CFT, Conf_nums, Conf_info,...
                        runopt, opt, args.num_replicas, stimulus, ...
                        'conf_DisplayName', conf_DisplayName, ...
                        'hfig', hfig2, ...
                        'plotflag', true, ...
                        'BinWidth', Time(1, 'ms'));
                end
            end
        end

        end
    end
end

%%

cmap = rwb( numel(FREQUENCY) );

[MM, mm] = deal(nan(numel(SIGNAL_TO_PLOT), numel(SR), ...
    numel(FREQUENCY), numel(AMPLITUDE)));

for r = 1:numel(SIGNAL_TO_PLOT)
    for s = 1:numel(SR)

        figure(111)
        hold on

        figure(222)
        hold on

        for j = 1:numel(FREQUENCY)
            
            for i = 1:numel(AMPLITUDE)
                unit = 's';
                Location = Time(RES.(SIGNAL_TO_PLOT{r}).(SR{s}).Location{i,j}, unit);
    
                f = Frequency(FREQUENCY(j), 'Hz');
                t = 1/f/2;
                t = Time(1, 'ms');
                
                [counts, edges] = histcounts( Location.ms, ...
                    'BinWidth', t.ms, ...
                    'Normalization', 'count');
                
                ff = 1/t;

                n_dyn_rep =  100;
    
                counts = ff.Hz * counts / args.num_replicas / n_dyn_rep;
    
                M = max(counts);
                m = mean(counts);
                MM(r,s,j,i) = M;
                mm(r,s,j,i) = m;
            end

            kw = {};
            kw = {'DisplayName', postprocess_window};
            if strcmp(postprocess_window, 'full')
                kw = [kw, {'Color', 'k'}];
            end
            % kw = {'Color', cmap(j,:)};

            figure(111)
            plot(AMPLITUDE, squeeze(MM(r,s,j,:)), kw{:});

            figure(222)
            plot(AMPLITUDE, squeeze(mm(r,s,j,:)), kw{:});
        end

        figure(111)
        title(sprintf('Max Rate %s %s', SR{s}, SIGNAL_TO_PLOT{r}));
        legend('Location', 'NorthWest')

        figure(222)
        title(sprintf('Mean Rate %s %s', SR{s}, SIGNAL_TO_PLOT{r}));
        legend('Location', 'NorthWest')

    end
end

%%

if do.inp_signal_phase

%     SIGNAL_TO_PLOT = {'Cleft'};
    % AMPLITUDE_TO_PLOT = cellfun(@(x) find(x==AMPLITUDE), {60});
    AMPLITUDE_TO_PLOT = 1:numel(AMPLITUDE);
    % FREQUENCY_TO_PLOT = cellfun(@(x) find(x==FREQUENCY), {200,2000,16000});
    FREQUENCY_TO_PLOT = 1:numel(FREQUENCY);
    SR_TO_PLOT = {'H1'};
    
    analysis.ant.inp_signal_phase(RES, variables, postprocess_window, FREQUENCY, AMPLITUDE, ...
        FREQUENCY_TO_PLOT, AMPLITUDE_TO_PLOT, SR_TO_PLOT, SIGNAL_TO_PLOT, runopt, opt)
    
end

%% number of data
if false
if exist('RES', 'var')
    for s = 1:numel(SR)
        NUM = NaN(numel(AMPLITUDE), numel(FREQUENCY));
        
        for i = 1:numel(AMPLITUDE)
            for j = 1:numel(FREQUENCY)
                if ~isempty(RES.(SIGNAL_TO_PLOT).(SR{s}).Location{i,j})
                    NUM(i,j) = sum(RES.Cleft.(SR{s}).Location{i,j} >= RES.analysisStartTime{i,j}.ms);
                end
            end
        end

        NUM
    end
end
end
end
