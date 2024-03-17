function [ resFiles, statistics ] = getMNAResult( stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt )
%GETMNARESULT Loads mechanical model results

if runopt.verbose >= 2
    fprintf('Loading MNA result...\n\t%s\n', runopt.path.oc_mna);
end

% try
    if runopt.recalculate.oc_mna == true
        overwrite = true;
    else
        overwrite = false;
    end
    
    t = tic();
    
    [resFiles, statistics] = oc_mna.ReorderResults( stimulus, topt, mechopt, mnaopt, runopt, opt, memopt, overwrite );
    
%     [ resFiles, runopt, stimulus, ~, mechopt, mnaopt, statistics ] = MNATransformResults( runopt.path.oc_mna, opt, memopt, overwrite );
% catch ME
%     warning(['Exception caught in reconstructResults, ', ... 
%              'setting `overwrite = true` and trying once again.']);
%     fprintf(['Exception details:\n', ...
%              '    identifier: %s\n', ...
%              '       message: %s\n'], ME.identifier, ME.message);
%     
%     overwrite = true;    
%     [ resFiles, runopt, stimulus, mechopt, mnaopt, statistics ] = MNATransformResults( runopt.path.oc_mna, opt, memopt, overwrite );
% end

    if strcmp(mechopt.middle_ear_params.identifier, 'none')
        resFiles.fMid_y = [];
        resFiles.fMid_t = [];
    else
        resFilesME = ResCacheMiddleEar( ...
            stimulus, topt, midopt, runopt, opt, memopt, paropt);
    
        resFiles.fMid_y = resFilesME.fy;
        resFiles.fMid_t = resFilesME.ft;
    end

    resFiles.options = matfile(fullfile(runopt.path.oc_mna, 'options.mat'), 'Writable', false);
    resFiles.settings = matfile(fullfile(runopt.path.oc_mna, 'settings.mat'), 'Writable', false);


if runopt.verbose >= 3

    load(fullfile(runopt.path.oc_mna, 'settings.mat'), ...
        'setStructure');

    % setStructure

end

dbprintf('... runtime %s.\n', disp_toc(toc(t)));

end

