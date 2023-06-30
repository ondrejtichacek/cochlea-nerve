function [ SynResFile ] = setSynapseResult( replicas, Voltage, Time, topt, antopt, runopt, paropt, memopt, SR, setStructure )
%SETSYNAPSERESULT Calls mechanical model, saves and returns its results

if ~exist(runopt.path.synapse, 'dir')
    mkdir(runopt.path.synapse);
end

if runopt.verbose >= 2
    fprintf('Synapse results will be saved to:\n%s\n', runopt.path.synapse);
end

save(fullfile(runopt.path.synapse, 'settings.mat'), ...
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


%% Handle result files (replicas)

fprintf('Simulating synapse with %d replicas ...\n', antopt.numberOfRepetitions);

SynResFile{length(replicas)} = [];

%%

% sr = SR{1};
sr = SR;

if paropt.useparalleltoolbox && ( paropt.useparfeval && ~paropt.submitAsBatchJob )

    tic

    V = parallel.pool.Constant(Voltage);
    T = parallel.pool.Constant(Time);
    
    for i = 1 : length(replicas)

        runopt_cpy = copy(runopt);
        runopt_cpy.verbose = min([2, runopt_cpy.verbose]);
        
        replicas(i).protect_PATH_error = true;
        
        f(i) = parfeval(p, @Synapse, 1, V, T, antopt, sr, replicas(i), topt, runopt_cpy, memopt);
    end
    fprintf('All jobs submited, %s\n', disp_toc);
    
else
    
    V = struct('Value', Voltage);
    T = struct('Value', Time);
    
end

if runopt.waitbarFunctionAvailable
    hw = waitbar(0, sprintf('Synapse: %d pending jobs ...\n', antopt.numberOfRepetitions*nSR), 'Resize', 'on');
end

tic
for i = 1 : length(replicas)
    
    if paropt.useparalleltoolbox && ( paropt.useparfeval && ~paropt.submitAsBatchJob )
        
        [idx,value] = fetchNext(f);
        value.protect_PATH_error = false;
        value.PATH = replicas(idx).dir; % transient property bug in parfeval workaround
        SynResFile{idx} = value;
        
    else
        replicas(i).protect_PATH_error = false;
        SynResFile{i} = Synapse(V, T, antopt, sr, replicas(i), topt, runopt, memopt);
        clear Synapse
    end
    
    progress_ratio = i/length(replicas);
    t_elapsed = toc;

    if runopt.waitbarFunctionAvailable
        waitbar( progress_ratio, hw, ...
            sprintf( ['Synapse: %d/%d/%d jobs completed\n', ...
                      'Elapsed time = %s\n', ...
                      'Est. remaining time = %s'], ...
                        i, ...
                        length(replicas), ...
                        antopt.numberOfRepetitions, ...
                        disp_toc(t_elapsed), ...
                        disp_toc(t_elapsed*(1/progress_ratio -1)) ) );
    end
    
    fprintf('%s Synapse sim. %d/%d/%d completed in %s, est. %s rem.\n', ...
        datestr(now), i, length(replicas), antopt.numberOfRepetitions, ...
        disp_toc(t_elapsed), disp_toc(t_elapsed*(1/progress_ratio -1)));
    
end

if runopt.waitbarFunctionAvailable
    close(hw);
end

%%

fprintf('Synapse simulation completed, %d reps, runtime: %s\n', length(replicas), disp_toc);

%% Create SUCCESS indicator

% fopen(fullfile(runopt.path.synapse, 'SUCCESS'),'w');

end
