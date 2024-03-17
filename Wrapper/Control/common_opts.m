function [args, opt, memopt, paropt] = common_opts(top_args, args)
arguments
    top_args % top-level args (either not defined here or overwriting the value)

    args.Numstacks (1,1) double = 600 

    args.computername (1,:) char = ''
    args.debug (1,1) logical = false
    args.gui (1,1) logical = false
    args.force_num_cores (1,1) double = 0
    
    args.memory (1,:) char = 'all'
    args.max_var_mem (1,:) char = ''
    
    args.noplot (1,1) logical = false
    args.parallel (1,1) logical = false
    args.parallel_switch = ''
    args.parallel_custom_monitor (1,1) logical = false
    args.parfeval_iteration_delay (1,1) Time = Time(0, 's')
    args.purge (1,1) logical = false
    args.recalculate (1,1) logical = false
    args.update_opts (1,:) cell = {}
    
    args.do_parts = []

    args.skip_do = {}
    args.do_recalculate = {}
    args.do_plot = {}
    args.do_override = struct()

    args.silent (1,1) logical = false
    args.do_not_change_settings (1,1) logical = false

    % ------------
    % stimulus

    args.product_args cell = {};

    args.Amplitude = []
    args.Frequency = []
    args.GlobalSamplingFrequency (1,1) Frequency = Frequency(200, 'kHz');
    args.amplitude_unit (1,:) char = 'spl' % 'phon'

    args.special_signal (1,:) char = ''
    args.signal_file (1,:) char = ''
    args.signal_var struct = struct();
    
    args.onset = 'quick' % char or Time
    args.offset = 'quick' % char or Time

    args.waveform (1,:) char = 'sine'

    args.zeroDuration (1,1) Time = Time(0.1, 'ms')
    args.fadeDuration (1,1) Time = Time(10, 'ms')

    args.t0 = Time(0);
    args.t2 = Time(0);
    args.tf_extra = 40; % 40 periods of the pure-tone
    args.strict_signal_length (1,1) logical = false
    
    args.stimulus_fun = []

    args.topt = []
    args.stimulus = []

    % ------------
    % outer and middle ear
    args.oe_identifier (1,:) char = 'Meddis'
    args.me_identifier (1,:) char = 'PBLL'


    % ------------

    args.mna_ver (1,:) char = 'dev2'
    args.ant_ver (1,:) char = 'v3'
    args.opt_ver (1,:) char = 'latest'
    
    % ------------
    % amplifier

    args.sim_groups cell = {'active', 'passive'};

    args.ohc_amplifier (1,:) char = 'electric'
    args.gain_factor (1,1) double = 1

    args.noise_damage (1,1) logical = false
    args.noise_damage_f0 (1,1) Frequency = Frequency(1200, 'Hz')
    
    % ------------
    % synapse

    args.ant_fun (1,:) char = 'ANT'

    args.slices = []

    args.SRS cell = {'H19b'}

    args.num_replicas (1,1) double = 1

    args.fiber_properties = struct([])

    args.calculate_only_ant_IC (1,1) logical = false

    % ------------
    % nerve

    args.hh_method (1,:) char = 'cw'

    % ------------

    args.analysis_compare_to = []
    args.analysis_plot_variables (1,:) cell = {}

    args.postprocess_windows_plot (1,:) cell = {}
    args.postprocess_window (1,:) char = ''
    
    args.time_limit (1,1) double = +Inf;  % in seconds
    args.ind_cx = 8;
    args.rng_seed (1,1) double = 0;
    
    args.no_create_only_load (1,1) logical = false%true
    args.no_create_only_load_oc_mna (1,1) logical = false%true
    args.cached_analysis (1,1) logical = false;
    args.cache = []
    
    args.heur_algorithm (1,:) char = 'DE'
    args.heur_algorithm_options cell = {}
    
    args.heur_res_src (1,:) char = ''
    args.heur_res_out (1,:) char = ''
    args.heur_res_cpt (1,:) char = ''
    
    args.heur_recalculate_population (1,1) logical = false
    args.parfeval_timeout = []

    args.population_size (1,1) int32 = 6
    args.maxiter (1,1) int32 = 1000
end

if ~isempty(top_args)
    args = mergestructs(top_args, args, 'nowarn', true);
end

if ~args.silent
    dispstruct(args)
end

if ~args.do_not_change_settings
    if args.debug == true
        % warning('on', 'all')
        dbstop if error
        
        % select workspace in the debug menu (Function Call Stack)
        % see variable by the openvar function
    else
        dbclear if error
    end
end
%%

rng(args.rng_seed);

%% Global options for the cochlear model

global computername

computername = args.computername;

if ~isempty(computername)
    opt = globalOpt([],computername);
else
    opt = globalOpt();
end

if args.force_num_cores ~= 0
    opt.num_cores = args.force_num_cores;
end

switch opt.computer
    case 'lucy'
        opt.resDir = '/qrk/CNaurum/res';
        opt.tmpDir = '/qrk/CNaurum/tmp';
        % opt.resDir = '/data/wrk/CNaurum/res';
        % opt.resDir = '/wrk/CNmagnesium/res';
        % opt.resDir = '/wrk/CNsodium/res';
    case 'cassandra'
        opt.resDir = '/data/CNaurum/res';
end

%% global arguments

% Memory options
% memopt = memOpt( ...
%     ...'numWorkers', 7, ... % if globalMemoryLimit (and numWorkers) is set, then maxVarMem & maxFileMem can be computed automaticaly
%     ...'globalMemoryLimit', memoryLimit, ...
%     ...'maxVarMem', 2000 * 1e6, ... % in Bytes
%     ...'maxFileMem', 2000 * 1e6, ...
%     ...'useAllAvailableMemory', true ... % true => maxVarMem & maxFileMem = current free memory
%     );

if ~isempty(args.max_var_mem)
    m = memOpt.parse_from_string(args.max_var_mem);
    memopt = memOpt( ...
        'useAllAvailableMemory', false, ...
        'maxVarMem', m.B, ...
        'maxFileMem', m.B); % in Bytes
else
    switch args.memory
        case 'all'        
            memopt = memOpt( ...
                'useAllAvailableMemory', true); % true => maxVarMem & maxFileMem = current free memory

        case 'auto'
            memopt = memOpt( ...
                'useAllAvailableMemory',false, ...
                'globalMemoryLimit', memoryLimit);

        otherwise

            m = memOpt.parse_from_string(args.memory);

            memopt = memOpt( ...
                'useAllAvailableMemory',false, ...
                'globalMemoryLimit', m.B); % in Bytes
    end
end    
        
% Parallel options
if strcmp(args.parallel_switch, 'heur') ...
    && ~strcmp(args.heur_algorithm, 'NR') ...
    && ~strcmp(args.heur_algorithm, 'TSEMO')
    
    paropt = parOpt( ...
        'useparalleltoolboxinheuristic', ...
            true && args.parallel, ...
            ...false, ...
        'useparalleltoolboxinouterscript', ...
            ...true && args.parallel, ...
            false, ...
        'useparalleltoolbox', ...
            ...true && args.parallel, ...
            false, ...
        'useparfor', false, ...
        'useparfeval', true, ...
        'submitAsBatchJob', false );
elseif strcmp(args.parallel_switch, 'one_job')
    paropt = parOpt( ...
        'useparalleltoolboxinouterscript', ...
            ...true && args.parallel, ...
            false, ...
        'useparalleltoolbox', ...
            ...false, ...
            true && args.parallel, ...
        'useparfor', false, ...
        'useparfeval', true, ...
        'submitAsBatchJob', false );
else
    paropt = parOpt( ...
        'useparalleltoolboxinouterscript', ...
            true && args.parallel, ...
            ...false, ...
        'useparalleltoolbox', ...
            ...true, ...
            false && args.parallel, ...
        'useparfor', false, ...
        'useparfeval', true, ...
        'submitAsBatchJob', false );
end

paropt.custom_monitor = args.parallel_custom_monitor;
paropt.force_num_cores = args.force_num_cores;

if args.parallel == true
   % paropt.configure_parcluster(paropt.force_num_cores);
   disp(paropt.parcluster);
end

end