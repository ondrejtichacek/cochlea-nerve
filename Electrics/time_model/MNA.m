function [ OC_results ] = MNA( ...
    stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, args )
%MNA
arguments
    stimulus 
    topt (1,1) timeOpt
    midopt (1,1) midOpt
    mechopt (1,1) mechOpt
    mnaopt (1,1) mnaOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt
    args.all_analysis_files_must_exist (1,1) logical = false
    args.variables_to_analyze = struct( ...
        'voltage', 'all', ...
        'current', 'all', ...
        'mech', 'all')
    args.analyse_signal_kwargs = struct();
end

t = tic;

% args.analyse_signal_kwargs = namedargs2cell(args.analyse_signal_kwargs);

if max([stimulus.frequency]) > 0
    mechopt.samplingFrequency = 64 * max([stimulus.frequency]);
    mnaopt.samplingFrequency = 64 * max([stimulus.frequency]);
end

%% Get MNA result

disp_title( 'Modified Nodal Analysis of Cochlear Network', '*' );

if runopt.ReconstructResults
    [ resFiles, runopt, stimulus, mechopt, mnaopt, statistics ] = MNAResult( ...
        stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt );
else
    MNAResult( stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt )
    variables = [];
    return
end

if runopt.no_create_only_load && ~runopt.found.oc_mna
    warning('LOAD_ONLY mode: MNA result not found! Setting variables to [] and trying to continue.');
    variables = [];

    args.all_analysis_files_must_exist = true;
end

%% Analysis

kwargs = namedargs2cell(args);
[variables, runopt] = oc_mna_analysis(resFiles, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, kwargs{:});

%%

OC_results = OCElectricsResult();
OC_results.result_files = resFiles;
OC_results.analysis = variables;


if runopt.delete_results_to_save_space
    % if args.nowarn == false
    %     warning('deleting results to save space - this should only be used in optimization etc.');
    % end

    fn = fieldnames(resFiles);

    fn = setdiff(fn, ...
        {'options', 'settings', ...
         'fMid_y', 'fMid_t'});

    for i = 1:numel(fn)
        f = resFiles.(fn{i});
        file_name = f.Properties.Source;
        if exist(file_name, 'file')
            f.Properties.Writable = true;
            vars = whos(f);
            for ii = 1:numel(vars)
                f.(vars(ii).name) = 'deleted';
            end
            % we need to free the space occupied by the variables
            mrepack(file_name);
        end
    end
end
%%

fprintf('MNA total runtime: %s\n', disp_toc(toc(t)))


end
