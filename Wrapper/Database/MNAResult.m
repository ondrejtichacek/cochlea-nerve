function [ resFiles, runopt, stimulus, mechopt, mnaopt, statistics ] = ...
    MNAResult( ...
        stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt )
%MNARESULT Loads or creates MNA model result
%   fMNAResult is a matfile with wariables
%       *
%       *
%       *
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
end

[runopt, ~] = findResult('mechanical', stimulus, [], midopt, mechopt, [], [], [], runopt, opt);

disp(runopt.hash.stimulus)

fprintf('Mechanical: %s\n', runopt.path.mech);

% [runopt, setStructure] = findResult('oc_mna', ...
%         stimulus, [], midopt, mechopt, mnaopt, [], [], runopt, opt, ...
%         'delete_on_error', false);

[runopt, setStructure] = findResultExtended(opt.get_extra_dirs(), 'oc_mna', ...
        stimulus, [], midopt, mechopt, mnaopt, [], [], runopt, opt, ...
        'delete_on_error', false);

fprintf('Electrical: %s\n', runopt.path.oc_mna);

if ~runopt.recalculate.oc_mna
    try
        topt_required = findTspan('electrical', runopt.found.oc_mna, runopt.path.oc_mna, topt, 'mnaopt', mnaopt);
    catch ME
        if strcmp(ME.message, 'Time points do not correspond.')
            warning('Time points do not correspond.')
            runopt.found.oc_mna = false;
            topt_required = findTspan('electrical', runopt.found.oc_mna, runopt.path.oc_mna, topt, 'mnaopt', mnaopt);
        else
            rethrow(ME)
        end
    end
end
topt.compute = copy(topt.total);

if runopt.no_create_only_load || runopt.no_create_only_load_oc_mna
    if ~runopt.found.oc_mna
        warning('LOAD_ONLY mode: MNA result not found!');
    end
    
elseif ... Calculate
        runopt.recalculate.oc_mna ...
        || ~runopt.found.oc_mna ...
        || ~isempty(topt_required)

    if runopt.recalculate.oc_mna
        action = 'new';
    else
        if runopt.found.oc_mna

            action = 'extend';
            topt.compute.tspan = topt_required.tspan;

            fprintf('Partial MNA result found:\n')
            fprintf('\t* tspan [%s, %s] \t ... found\n', topt.total.t0, topt.compute.t0)
            fprintf('\t* tspan [%s, %s] \t ... will be computed\n', topt.compute.t0, topt.compute.tf)

        else
            action = 'new';
        end
    end

    setMNAResult( action, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, setStructure );

end

if runopt.no_create_only_load && ~runopt.found.oc_mna && ~runopt.no_create_only_load_oc_mna
    warning('LOAD_ONLY mode: MNA result not found!');
    [resFiles, statistics ] = deal([]);
    
elseif nargout > 0 || runopt.ReconstructResults == true

    if nargout > 0 && runopt.ReconstructResults == false

        warning('MNA results must be reconstructed first!');

    end

    [ resFiles, statistics ] = getMNAResult( stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt);
%     [resFiles, runopt_old, stimulus, mechopt, mnaopt, statistics] = getMNAResult( runopt, opt, memopt );

%     runopt_old.updateWithObject( runopt ); % update values in runopt_old by values in runopt
%     runopt = runopt_old;

end

end
