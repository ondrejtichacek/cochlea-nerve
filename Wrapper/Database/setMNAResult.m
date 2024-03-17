function setMNAResult( ...
    action, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, setStructure )
%SETMNARESULT Calls electrical model, saves and returns its results

%%

tic_MNA = tic;

if ~exist(runopt.path.oc_mna, 'dir')
    mkdir(runopt.path.oc_mna);
end

%% Check for lock file

if runopt.clear_lock_files_if_found == true
    handleLockFile(runopt.path.oc_mna, 'action', 'purge');
end

handleLockFile(runopt.path.oc_mna, 'action', 'check');

%% Cleanup -- delete old result files if exist

switch action
    case 'new'

        deleteResultFiles(runopt.path.oc_mna);

        y0 = [];

    case 'extend'

        cp_t = load(fullfile(runopt.path.oc_mna, 'cp_t.mat'), 't', 'unit');

        i = find(Time(cp_t.t, cp_t.unit) == topt.compute.t0, 1);
        assert(~isempty(i), 'Requested time step %s not found in checkpoints.mat.', topt.compute.t0);
        
        cp_y = matfile(fullfile(runopt.path.oc_mna, 'cp_Y.mat'), 'Writable', false);
        
        y0 = cp_y.y(i,:)';

        clear('cp_t', 'cp_y')

    otherwise
        error('Unknown option %s.', action);
end

%%

lock_file_cleanup = handleLockFile(runopt.path.oc_mna, 'action', 'create');

runopt = oc_mna.mainMNA( stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, y0 );


%%

save( fullfile( runopt.path.oc_mna, 'options.mat'), ...
    'stimulus', ...
    'topt', ...
    'midopt', ...
    'mechopt', ...
    'mnaopt', ...
    'runopt', ...
    'opt', ...
    'memopt', ...
    'paropt', ...
    '-v7.3');


save(fullfile(runopt.path.oc_mna, 'settings.mat'), ...
    'stimulus', 'topt', 'mechopt', 'memopt', 'setStructure', ...
    '-v7.3');

fprintf('\nThe MNA model part evaluated in %s\n', disp_toc(toc(tic_MNA)));


%% Create SUCCESS indicator

fid = fopen(fullfile(runopt.path.oc_mna, 'SUCCESS'),'w');
fclose(fid);

clear lock_file_cleanup

end

