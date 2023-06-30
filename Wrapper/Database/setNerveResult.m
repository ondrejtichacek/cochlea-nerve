function [ NerveResFile, runopt ] = setNerveResult( replicas, stimulus, topt, mechopt, mnaopt, antopt, hhopt, opt, runopt, memopt, paropt, setStructure  )
%SETNERVERESULT Calls nerve model, saves and returns its results

tic

mkdir(runopt.path.nerve);

if runopt.verbose >= 2
    fprintf('Nerve results will be saved to:\n%s\n', runopt.path.nerve);
end

save(fullfile(runopt.path.nerve, 'settings.mat'), ...
    'setStructure', 'topt');



if paropt.useparalleltoolbox && ( paropt.useparfeval && ~paropt.submitAsBatchJob )
    tic
    % p = gcp();
    p = paropt.startparpool();
    fprintf('Parpool created (loaded) in %s.\n', disp_toc);
    
    % shuffle local RNG, seed RNG on workers, set worker RNG algorithm to Mersene-Twister
    rngParallel('twister','shuffle',p);
else
    rng('shuffle','twister')
end


%%

slices = 1:length(antopt.slicesToExcite);

%% Handle result files (replicas)

fprintf('Simulating nerve with %d replicas ...\n', antopt.numberOfRepetitions);


numRepToDo = numel(replicas);

if runopt.waitbarFunctionAvailable
    hw = waitbar(0, sprintf('Nerve: %d pending jobs ...\n',antopt.numberOfRepetitions),'Resize','on');
end

%%

if paropt.useparalleltoolbox && ( paropt.useparfeval && ~paropt.submitAsBatchJob )
    
    tic
    
    for i = 1:numRepToDo

        runopt_cpy = copy(runopt);
        runopt_cpy.verbose = min([2, runopt_cpy.verbose]);
        
        replicas(i).protect_PATH_error = true;
        replicas(i).src_data.protect_PATH_error = true;
        
        scaleFactor = hhopt.scaleFactor;
        f(i) = parfeval(p, @cochleaRunBruce, 1, ...
                        topt, ...
                        hhopt, ...
                        scaleFactor, slices, ...
                        replicas(i), runopt_cpy, memopt);

    end

    fprintf('All jobs submited, %s\n', disp_toc);

end


%%

NerveResFile{numRepToDo} = [];
for i = 1:numRepToDo
    
    
    if paropt.useparalleltoolbox && ( paropt.useparfeval && ~paropt.submitAsBatchJob )

        timeout = 900; % DEBUG
        timeout = 0;
        if timeout > 0
            [idx, value] = fetchNext(f, timeout);
        else
            [idx, value] = fetchNext(f);
        end
        value.PATH = replicas(idx).dir; % transient property bug in parfeval workaround
        value.protect_PATH_error = false;
        
    else
        
        idx = i;
        
        replicas(idx).protect_PATH_error = false;
                
        scaleFactor = hhopt.scaleFactor;
        value = cochleaRunBruce( ...
                    topt, ...
                    hhopt, ...
                    scaleFactor, slices, ...
                    replicas(idx), runopt, memopt);
        
    end
        

    NerveResFile{idx} = value;    
    
    progress_ratio = i / antopt.numberOfRepetitions;
    t_elapsed = toc;

    if runopt.waitbarFunctionAvailable
        waitbar( progress_ratio, hw, ...
            sprintf( ['Nerve: %d/%d jobs completed\n',...
                      'Elapsed time = %s\n', ...
                      'Est. remaining time = %s'], ...
                        i, ...
                        antopt.numberOfRepetitions, ...
                        disp_toc(t_elapsed), ...
                        disp_toc(t_elapsed*(1/progress_ratio -1)) ) );
    end
    
    fprintf('%s nerve simulation %d/%d complete.\n', datestr(now), i, antopt.numberOfRepetitions);
    
end

if runopt.waitbarFunctionAvailable
    close(hw)
end

texec = toc;

fprintf('Nerve simulation completed, %d reps, runtime: %s\n', antopt.numberOfRepetitions, disp_toc(texec));

%% runopt

runopt.update('texec_nerve', texec);

%% Create SUCCESS indicator

% fopen(fullfile(runopt.path.nerve, 'SUCCESS'),'w');

end
