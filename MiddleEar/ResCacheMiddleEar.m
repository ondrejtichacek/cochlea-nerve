function [ resFiles ] = ResCacheMiddleEar( ...
    stimulus, topt, midopt, runopt, opt, memopt, paropt)
%RESCACHEMIDDLEEAR
arguments
    stimulus (1,:) {mustBeA(stimulus, ["Signal", "CompoundSignal"])}
    topt     (1,1) timeOpt
    midopt   (1,1) midOpt
    runopt   (1,1) runOpt
    opt      (1,1) globalOpt
    memopt   (1,1) memOpt
    paropt   (1,1) parOpt
end

%% Get Middle Ear result

disp_title( 'Middle Ear - Modified Nodal Analysis', '*' );

[runopt, setStructure] = findResult('middleear', stimulus, [], midopt, [], [], [], [], runopt, opt);

fprintf('Middle Ear: %s\n', runopt.path.mid);

try
    topt_required = findTspan('middleear', runopt.found.mid, runopt.path.mid, topt, 'midopt', midopt);
catch ME
    if strcmp(ME.message, 'Time points do not correspond.')
        warning('Time points do not correspond.')
        runopt.found.mid = false;
        topt_required = findTspan('middlear', runopt.found.mid, runopt.path.mid, topt, 'midopt', midopt);
    else
        rethrow(ME)
    end
end
topt.compute = copy(topt.total);

if runopt.no_create_only_load && ~runopt.found.mid
    warning('LOAD_ONLY mode: Middle Ear result not found!');
    [ft, fy] = deal([]);
    
elseif any([... %Calculate
        ~runopt.found.mid, ...
        runopt.recalculate.mid, ...
        ~isempty(topt_required)])

    if runopt.recalculate.mid
        action = 'new';
    else
        if runopt.found.mid

            action = 'extend';
            topt.compute.tspan = topt_required.tspan;

            fprintf('Partial Middle Ear result found:\n')
            fprintf('\t* tspan [%s, %s] \t ... found\n', topt.total.t0, topt.compute.t0)
            fprintf('\t* tspan [%s, %s] \t ... will be computed\n', topt.compute.t0, topt.compute.tf)

        else
            action = 'new';
        end
    end

    [ft, fy] = setMiddleEearMNAResult(action, stimulus, topt, midopt, runopt, opt, memopt, paropt, setStructure );

elseif runopt.found.mid
    
    [ft, fy] = getMiddleEearMNAResult(runopt);
    
end

resFiles = struct( ...
    'ft', ft, ...
    'fy', fy);

end

function [ft, fy] = setMiddleEearMNAResult( ...
    action, stimulus, topt, midopt, runopt, opt, memopt, paropt, setStructure )
arguments
    action   (1,:) char
    stimulus (1,1) {mustBeA(stimulus, ["Signal", "CompoundSignal"])}
    topt     (1,1) timeOpt
    midopt   (1,1) midOpt
    runopt   (1,1) runOpt
    opt      (1,1) globalOpt
    memopt   (1,1) memOpt
    paropt   (1,1) parOpt
    setStructure   struct
end
%% Cleanup -- delete old result files if exist

switch action
    case 'new'

        deleteResultFiles(runopt.path.mid);

        y0 = [];

    case 'extend'

%         cp_t = load(fullfile(runopt.path.mid, 'cp_t.mat'), 't', 'unit');
% 
%         i = find(Time(cp_t.t, cp_t.unit) == topt.compute.t0, 1);
%         assert(~isempty(i), 'Requested time step %s not found in checkpoints.mat.', topt.compute.t0);
%         
%         cp_y = matfile(fullfile(runopt.path.mid, 'cp_Y.mat'), 'Writable', false);
%         
%         y0 = cp_y.y(i,:)';

        settings = load(fullfile(runopt.path.mid, 'settings.mat'), 'midopt');
        unit_t = settings.midopt.simulation_units.time;

        ft = load(fullfile(runopt.path.mid, 't.mat'), 't');

        tol = Time(1e-9, 'ms');
        i = find(abs(Time(ft.t, unit_t) - topt.compute.t0) < tol, 1);
        assert(~isempty(i), 'Requested time step %s not found in checkpoints.', topt.compute.t0);
        
        fy = matfile(fullfile(runopt.path.mid, 'y.mat'), 'Writable', false);
        
        y0 = fy.y(i,:)';

        
        clear('ft', 'fy')

    otherwise
        error('Unknown option %s.', action);
end

%%

tic_Middle_Ear = tic;

if ~exist(runopt.path.mid, 'dir')
    mkdir(runopt.path.mid);
end

%% Evaluate model

[t, y, y_desc] = MiddleEar_MNA(topt, stimulus, midopt, runopt, opt);

%% Save results

if runopt.do_not_save_results == true
    % Sometimes we don't want to save results to disk (e.g. in heuristic
    % optimization). Keep them in memory for now
    
    ft = struct('t', t);
    fy = struct('y', y, 'desc', y_desc);
    
else
    save( fullfile(runopt.path.mid, 'options.mat'), ...
        'stimulus', ...
        'topt', ...
        'midopt', ...
        'runopt', ...
        'opt', ...
        'memopt', ...
        'paropt', ...
        '-v7.3');

    save(fullfile(runopt.path.mid, 'settings.mat'), ...
        'stimulus', 'topt', 'midopt', 'memopt', 'setStructure', ...
        '-v7.3');

    ft = matfile(fullfile(runopt.path.mid, 't.mat'), 'Writable', true);
    fy = matfile(fullfile(runopt.path.mid, 'y.mat'), 'Writable', true);

    ft.t = t;
    fy.y = y;
    fy.desc = y_desc;

    % Create SUCCESS indicator
    fid = fopen(fullfile(runopt.path.mid, 'SUCCESS'),'w');
    fclose(fid);
end

fprintf('\nThe Middle Ear model part evaluated in %s\n', disp_toc(toc(tic_Middle_Ear)));

end

function [ ft, fy ] = getMiddleEearMNAResult(runopt)
arguments
    runopt (1,1) runOpt
end

if runopt.verbose >= 2
    fprintf('Loading Middle Ear result...\n\t%s\n', runopt.path.mid);
end

ft = matfile(fullfile(runopt.path.mid, 't.mat'), 'Writable', true);
fy = matfile(fullfile(runopt.path.mid, 'y.mat'), 'Writable', true);

if runopt.verbose >= 3

    load(fullfile(runopt.path.mid, 'settings.mat'), ...
        'setStructure');

    display(setStructure)

end

end




