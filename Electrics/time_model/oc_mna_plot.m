function [] = oc_mna_plot( plotopt, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt, results)
arguments
    plotopt (1,1) plotOpt
    stimulus 
    topt (1,1) timeOpt
    midopt (1,1) midOpt
    mechopt (1,1) mechOpt
    mnaopt (1,1) mnaOpt
    antopt
    hhopt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt
    results
end
%OC_MNA_PLOT 


if runopt.plot.oc_mna == true ...
        && plotopt.doplot ...
        && ( ...
            plotopt.interactive_selection ...
            || any(structfun(@(x) x, plotopt.do.mechanical)) ...
            || any(structfun(@(x) x, plotopt.do.oc_electric)))
    
    [ resFiles, statistics ] = oc_mna.ReorderResults(stimulus, topt, mechopt, mnaopt, runopt, opt, memopt);
    
    results.oc_mna = resFiles;

    if strcmp(mechopt.middle_ear_params.identifier, 'none')
        resFilesME = [];
    else
        resFilesME = ResCacheMiddleEar( ...
            stimulus, topt, midopt, runopt, opt, memopt, paropt);
    end

    % resFiles.fMid_y = resFilesME.fy;
    % resFiles.fMid_t = resFilesME.ft;

    results.middle_ear = resFilesME;

    if ~isfield(results, 'synapse')
        results.synapse = [];
    end
    if ~isfield(results, 'nerve')
        results.nerve = [];
    end

    if plotopt.interactive_selection == true
        
        plotfcn = @(plotopt)  ...
            MEplotting(results, plotopt, stimulus, topt, midopt, runopt, opt, memopt) ...
            & ...
            MechPlotting(results, plotopt, stimulus, topt, midopt, mechopt, runopt, opt, memopt) ...
            & ...
            MNAplotting(results, plotopt, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt);
        
        plotopt.InteractiveSelection(plotfcn, 'middle_ear', 'mechanical', 'oc_electric');
    else
        MechPlotting(results, plotopt, stimulus, topt, midopt, mechopt, runopt, opt, memopt);
        MNAplotting(results, plotopt, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt);
    end
    
% elseif runopt.ReconstructResults == true
% 
%     MNATransformResults(runopt.path.oc_mna, opt, memopt);

end

end

