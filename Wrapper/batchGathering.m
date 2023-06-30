function [ statistics ] = ...
        batchGathering(ANALYSIS_FUN, n_arg_out, Configurations, ...
            stimulus, SIM_OPTS, ...
            runopt, opt, memopt, paropt, ...
            args)

arguments
    ANALYSIS_FUN
    n_arg_out (1,1) %int32
    Configurations (1,:) struct
    
    stimulus (1,1) function_handle
    SIM_OPTS
    
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt
    
    args.override_analysis_start = []
    args.verbose (1,1) int32 = 2
    args.debug (1,1) logical = false
    args.DRY_RUN (1,1) logical = false
    args.fig_save_dir_name_fun = @(Conf) fullfile(opt.cochleadir, 'Results', sprintf( ...
            '%dHz_%ddB', Conf.Frequency, Conf.Amplitude));
    args.PLOT_FUN = [];
    args.plotopt = [];
    
    args.pass_extended_options (1,1) logical = true
    
    args.analysis_fun_kwargs = {};
    
    args.skip_configuration_uniqueness_test (1,1) logical = false
    
    args.interactive_stop (1,1) logical = ~(usejava('jvm') && ~feature('ShowFigureWindows'))

    args.parfeval_iteration_delay (1,1) Time = Time(0, 's')
end

t.analysis = tic;

if args.verbose
    print = @(varargin) fprintf(varargin{:});
else
    print = @(varargin) disp('');
end

%% Check configurations

print("The configuration list contains the following variables:\n")
names = fieldnames(Configurations);
for i = 1:numel(names)
    name = names{i};
    switch class([Configurations.(name)])
        case 'double'
            m = min([Configurations.(name)]);
            M = max([Configurations.(name)]);
            if isa(m, 'Unit')
                print("    %s: range [%s, %s]\n", name, m, M)
            else
                print("    %s: range [%g, %g]\n", name, m, M)
            end
        case 'cell'
            print("    %s: values:", name)
            disp("for now not printing")
            % disp(unique([Configurations.(name)]))
    end
end

print("Checking if configurations are unique ...\n");
hashes = cell(numel(Configurations),1);
for i = 1:numel(Configurations)
    hashes{i} = DataHash(Configurations(i));
end
uhashes = unique(hashes);
if numel(uhashes) ~= numel(hashes)
    error("Configurations are not unique!")
end

print("We will run %d simulations.\n", numel(Configurations))

%% Create parallel pool if required

if paropt.useparalleltoolboxinouterscript
     p = paropt.startparpool(opt.num_cores);
end

%% Submit jobs (and gather data)

statistics = cell(numel(Configurations), n_arg_out);
hashes = cell(numel(Configurations), 1);

if args.skip_configuration_uniqueness_test
    warning("Skipping configuration uniqueness test #2")
end

resFiles = {};

% conf_parameters = fieldnames(Configurations);

for i = 1:numel(Configurations)
    
    Conf = Configurations(i);
    
    [runopt_cpy, opt_cpy, memopt_cpy, paropt_cpy] = copy_all( ...
            runopt, opt, memopt, paropt);

    SIM_OPTS_cpy = cell(size(SIM_OPTS));
    for j = 1:numel(SIM_OPTS)
        if isa(SIM_OPTS{j}, 'Opt')
            SIM_OPTS_cpy{j} = copy(SIM_OPTS{j});
        elseif isa(SIM_OPTS{j}, 'function_handle')
            SIM_OPTS_cpy{j} = SIM_OPTS{j}(Conf);
        else
            error('Expected the %d argument of SIM_OPTS to be an instance of Opt (or children classes) or a function handle, not a %d', j, class(SIM_OPTS{j}))
        end
    end
    
    % conf_par = cell(numel(conf_parameters),1);
    % for ii = 1:numel(conf_parameters)
    %     par = conf_parameters{ii};
    %     conf_par{ii} = Conf.(par);
    % end

    % create time opt & stimulus
    % [topt_cpy, stimulus_cpy] = stimulus(conf_par{:});
    [topt_cpy, stimulus_cpy] = stimulus(Conf);

    zero_time = stimulus_cpy.zeroDuration;
    analysis_start = stimulus_cpy.default_analysis_start_time;

    if ~isempty(args.override_analysis_start)
        print('Overriding analyss start time from %s to %s\n', analysis_start, args.override_analysis_start);
        analysis_start = args.override_analysis_start;
    end
    
    runopt_cpy.update( ...
        'zeroTime', zero_time, ...
        'analysisStart', analysis_start, ...
        'figureSaveDir', args.fig_save_dir_name_fun(Conf));

    if paropt.useparalleltoolboxinouterscript
        paropt_cpy.useparalleltoolbox = false;
    end        
    
    if args.pass_extended_options
        OPT = [{stimulus_cpy, topt_cpy}, ...
                SIM_OPTS_cpy, ...
               {runopt_cpy, opt_cpy, memopt_cpy, paropt_cpy}];
    else
        OPT = [{stimulus_cpy, topt_cpy}, ...
                SIM_OPTS_cpy];
    end
    
    % check uniqueness again
    if args.skip_configuration_uniqueness_test
    else
        % hashes{i} = DataHash([{stimulus_cpy, topt_cpy}, SIM_OPTS_cpy]);
        hashes{i} = DataHash([{stimulus_cpy}, SIM_OPTS_cpy]);
        if i > 1
            if ~isempty(find(strcmp(hashes{i}, hashes(1:i-1)),1))
                error("This simulation configuration does not seem to be unique")
            end
        end
    end
     
    if args.DRY_RUN
        printf('Dry run %d:', i)
        disp(topt_cpy.total)
    else
        if paropt.useparalleltoolboxinouterscript
            pf(i) = parfeval(p, ANALYSIS_FUN, n_arg_out, OPT{:}, args.analysis_fun_kwargs{:});

            if i < numel(Configurations)
                if args.parfeval_iteration_delay.s > 0
                    dbprintf('Pausing %s before submitting the next iteration.', args.parfeval_iteration_delay)
                    pause(args.parfeval_iteration_delay.s)
                end
            end
        else
%             try
                disp(topt_cpy.total)
                [statistics{i,:}] = ANALYSIS_FUN(OPT{:}, args.analysis_fun_kwargs{:});
%             catch ME
%                 if any(strcmp(ME.identifier, ...
%                         {'MATLAB:minrhs', ... % 'Not enough input arguments.'
%                          'MATLAB:dispatcher:InexactCaseMatch', ...
%                         }))
%                     rethrow(ME)
%                 end
%                 
%                 warning('Something went wrong with the job #%d',i)    
%                 disp(getReport(ME, 'extended'))
% 
%                 warning('Deleting folders\n%s\n%s\ncontents and trying again', ...
%                     runopt_cpy.path.synapse, runopt_cpy.path.nerve)                        
%                 srmdir(runopt_cpy.path.synapse, 's')
%                 srmdir(runopt_cpy.path.nerve, 's')
% 
%                 statistics{i} = ANALYSIS_FUN(OPT{:});
%             end
            print('Iteration %d/%d done \n', i, numel(Configurations))
            
            if ~isempty(args.plotopt) && ~isempty(args.PLOT_FUN)
                if args.plotopt.doplot == true
                    args.PLOT_FUN(args.plotopt, OPT{:});
                end
            end
        end
    end 
end

if args.DRY_RUN
    warning('Ending the dry run')
    return
end

%% Gather data (if parallel)

if paropt.useparalleltoolboxinouterscript
    if args.verbose == 1
        fprintf('\n %d total:   [', numel(Configurations))
    end
    
    if paropt.custom_monitor == true
        fetchNextMonitor(pf);
    end
    
    if args.interactive_stop == true
%         Hf = figure;
%         H = uicontrol( ...
%             'Style', 'PushButton', ...
%             'String', 'Break batchGathering - fetchNext', ...
%             'Callback', 'delete(gcbf)');

        pfh = parfeval_ui_monitor(pf, Configurations);

        cleanup = onCleanup(@()myCleanupFun(pfh.figure));
        
        fetch_next_timeout = {60};
    else
        fetch_next_timeout = {};
    end
    
    for i = 1:numel(Configurations)
        try
            wait_for_job = true;
            
            while wait_for_job
                value = cell(1,n_arg_out);
                [idx, value{1,:}] = fetchNext(pf, fetch_next_timeout{:});
                [statistics{idx,:}] = value{1,:};
                
                if ~isempty(idx) % successful fetch
                    wait_for_job = false;
                else
                    if args.interactive_stop == true
                        pfh = parfeval_ui_monitor(pf, Configurations, 'handles', pfh);
                        H = pfh.stop_button;
                        if not(ishandle(H))
                            % stop_reason = 'user';
                            % do_stop = true;
                            delete(gcp);
                            return
                        end
                    end
                end
            end
        catch ME
            warning('Something went wrong with the job #%d', i)
                disp(Configurations(i))
                disp(pf(i))
                disp(getReport(ME, 'extended'))
            continue
        end
        if args.verbose == 1
            fprintf(' #%d (%s),', i, disp_toc(toc(t.analysis)))
        elseif args.verbose == 2
            fprintf('\n    %d/%d completed, %s,    %s', ...
                i, numel(Configurations), disp_toc(toc(t.analysis)), ...
                pf(idx).InputArguments{1}.name)
        end
    end
    if args.verbose == 1
        fprintf(' ]\n')
    end
end

if args.verbose
    fprintf('\nAnalysis (loading) completed in %s\n', disp_toc(toc(t.analysis)));
end

    function myCleanupFun(Hf)
        if ~isempty(Hf)
            close(Hf);
        end
    end


end