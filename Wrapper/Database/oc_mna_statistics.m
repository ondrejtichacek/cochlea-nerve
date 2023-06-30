function [ OC_results ] = oc_mna_statistics( ...
        stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, args )
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
    args.all_analysis_files_must_exist (1,1) logical = true
    args.variables_to_analyze = struct( ...
        'voltage', 'all', ...
        'current', 'all', ...
        'mech', 'all')
    args.analyse_signal_kwargs = struct();
end

if max([stimulus.frequency]) > 0
    mechopt.samplingFrequency = 64 * max([stimulus.frequency]);
    mnaopt.samplingFrequency = 64 * max([stimulus.frequency]);
end

[runopt] = findResultExtended(opt.get_extra_dirs(), 'oc_mna_analysis', ...
        stimulus, [], midopt, mechopt, mnaopt, [], [], runopt, opt, ...
        'delete_on_error', false);

resFiles = [];
kwargs = namedargs2cell(args);
[variables, runopt] = oc_mna_analysis(resFiles, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, kwargs{:});

OC_results = OCElectricsResult();
% OC_results.result_files = resFiles;
OC_results.analysis = variables;


end
