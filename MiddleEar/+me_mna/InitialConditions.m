function [Y0_consistent, Yp0_consistent] = InitialConditions( ...
    ODEFUN, solveropt, midopt, runopt, opt, y0, Y0, args)
%INITIALCONDITIONS
arguments
    ODEFUN
    solveropt
    midopt 
    runopt 
    opt 
    y0 
    Y0
    args.max_tries = 100
end

solveropt_copy = odeset(solveropt, ...
    'RelTol', 1e-9);

n = numel(Y0);

size_info = struct(...
    'mid', struct( ...
        'start', 1, ...
        'end', n, ...
        'size', n ...
    ));

extra_dirs = {runopt.read_only_cache};

[runopt, ~] = findResultExtended(extra_dirs, 'me_mna_dae_ic', ...
        [], [], midopt, [], [], [], [], runopt, opt);

cache_dae_file = fullfile(runopt.path.me_mna_dae_ic, 'dae_ic.mat');
if ~exist(runopt.path.me_mna_dae_ic, 'dir')
    mkdir(runopt.path.me_mna_dae_ic)
end

% Compute consistent inital conditions of the DAE for the solver
% if not continuing an existing simulation
if ~runopt.recalculate.me_mna_dae_ic ...
        && isempty(y0) ...
        && runopt.found.me_mna_dae_ic ...
        && exist(cache_dae_file, 'file')
    
    dbprintf('Loading consistent initial conditions from a file:\n -- %s.\n', cache_dae_file);
    tmp = load(cache_dae_file);
    y0 = tmp.y0;
    yp0 = tmp.yp0;
    
elseif runopt.no_create_only_load_me_mna_dae_ic
    
    error('We require initial conditions to exist.')
    
elseif isempty(y0)
    
    OutputFcnList = {};

    dbprintf('Computing consistent initial conditions.\n');
    timerInitialConditions = tic;
    
    try_again = true;
    try_num = 0;
    ic_system_method = 'M_full';
    
    while try_again
        try_num = try_num + 1;
        
        dbprintf(' -- Using ''%s'' modification of the system to compute initial conditions\n', ic_system_method)
        
        switch ic_system_method
            case 'M_full'
                solveropt_init = solveropt_copy;
            otherwise
                rethrow(E)
        end

        if try_num == 1
            tspan = Time([0, 10], 'us'); % super short tspan allows IC to pass daeic3
        elseif try_num < 10
            tspan = Time([0, 1000], 'ms');
        else
            tspan = Time([0, 1000], 'ms');
        end

        dbprintf(' -- tspan: [%s, %s]\n', tspan(1), tspan(2))

        tsp = tspan.s;
        
        OutputFcnList_init = OutputFcnList;

        if runopt.verbose >= 3
        % Display solver statistics
        % Waitbar function as OutputFcn, to be used after every integration timestep
        
            ws = WaitbarStorage();
            ws.num_div = 1;
        
            if runopt.waitbarFunctionAvailable
                OutputFcnList_init{end+1} = @(t, y, flag) odewaitbar(t, y, flag, tsp(end), ws);
            else
                OutputFcnList_init{end+1} = @(t, y, flag) odewaittext(t, y, flag, tsp(end), ws);
            end
        end
        
        if numel(OutputFcnList_init) > 0
            OutputFcn = @(t,y,flag) OutputFcnChain(OutputFcnList_init, t, y, flag);
            
            % number of ~ must correspond with the number of extra arguments in ode15s
            OutputFcnWrapper = @(t, y, flag, ~, ~, ~, ~) OutputFcn(t, y, flag);
            solveropt_init = odeset(solveropt_init, 'OutputFcn', OutputFcnWrapper);
        end
        
        Y0_tmp = Y0;
            
        if try_num > 1

            if false % try_num >= 10

                if ~isempty(solveropt_init.Jacobian)
                    solveropt_init = odeset(solveropt_init, ...
                        'Jacobian', JACOBFUN);
                end
                
                mechopt.vohc_0 = vars.vohc;

                m = size_info.mech.start;
                n = size_info.mech.end;
                Y0_tmp(m:n) = 0;
                yp0(m:n) = 0;
            end
            
            solveropt_init = odeset(solveropt_init, 'InitialSlope', yp0);
            solveropt_init = mtode.odeset(solveropt_init, 'SkipICCheck', true);
        end
        
        timerInitialConditions_iter = tic;

        % try
        sol = mtode.ode15s(ODEFUN, tsp, Y0_tmp, ...
                mtode.odeset(solveropt_init,  ...
                    'Stats', 'off'));
        % catch E
        %     if strcmp(E.identifier, 'MATLAB:daeic3:NeedBetterY0')
        %         dbprintf('Computing consistent IC failed with method %d\n', try_num)
        %         continue
        %     else
        %         rethrow(E)
        %     end
        % end

        % yy = sol.y';            
        % y0 = yy(end,:)';
        % yp0 = sol.idata.dif3d(:,end,end);
        
        [y0, yp0] = mtode.deval(sol, tsp(end));
        
        dbprintf('initial condition increment norm diff:\n')
        step_norm.total = norm(Y0 - y0);
        norm_change.total = norm(Y0 - y0)/norm(Y0);
        Z = Y0 - y0;
        Z(abs(Z) < 1e-12) = 0;
        adj_norm_change.total = norm(Z)/norm(Y0);
        dbprintf(' -- -- %s   %g   %g %%   %g %%\n', 'total', step_norm.total, 100*norm_change.total, 100*adj_norm_change.total)
        fn = fieldnames(size_info);
        for k = 1:numel(fn)
            
            if strcmp(fn{k}, 'blocks')
                continue
            end
            
            m = size_info.(fn{k}).start;
            n = size_info.(fn{k}).end;
            
            NN = norm(Y0(m:n));
            if NN == 0
                step_norm.(fn{k}) = 0;
                norm_change.(fn{k}) = 0;
                adj_norm_change.(fn{k}) = 0;
            else
                step_norm.(fn{k}) = norm(Y0(m:n) - y0(m:n));
                norm_change.(fn{k}) = norm(Y0(m:n) - y0(m:n))/NN;
                Z = Y0(m:n) - y0(m:n);
                Z(abs(Z) < 1e-12) = 0;
                adj_norm_change.(fn{k}) = norm(Z)/NN;
            end
            
            dbprintf(' -- -- %s   %g   %g %%   %g %%\n', fn{k}, step_norm.(fn{k}), 100*norm_change.(fn{k}), 100*adj_norm_change.(fn{k}))

        end
        Y0 = y0;

        dbprintf(' -- iteration runtime %s.\n', disp_toc(toc(timerInitialConditions_iter)));
        
        if try_num <= 20
            try_again = true;
        elseif try_num >= args.max_tries
            try_again = false;
            warning('Maximum number of iterations exceeded, finishing initial conditions optimization.')
        elseif try_num > 20
            if all(structfun(@(x) x < 1e-6, norm_change)) || ...
               all(structfun(@(x) x < 1e-9, adj_norm_change))
                try_again = false;
            end
        end
        
        ic_system_method = 'M_full';
        N = numel(tspan);
    end    
    
    %% Sanity checking
    
    % max_expected_values = struct( ...
    %     'mid', 1e-9);
    % 
    % fn = fieldnames(max_expected_values);
    % for k = 1:numel(fn)
    % 
    %     m = size_info.(fn{k}).start;
    %     n = size_info.(fn{k}).end;
    % 
    %     if max(abs(Y0(m:n))) > max_expected_values.(fn{k})
    %         error('Initial conditions did not converge! Maximal expected value of %g for %s exceeded. Value = %g', ...
    %             max_expected_values.(fn{k}), fn{k}, max(abs(Y0(m:n))))
    %     end
    % end
    
    %% Finalizing
    
    dbprintf('... runtime %s.\n', disp_toc(toc(timerInitialConditions)));
    % figure; surf(getVoltage(yy, 'IHC', mnaopt))
    
    fprintf('Saving initial guess for MNA DAE to %s\n', cache_dae_file)
    save(cache_dae_file, 'y0', 'yp0', '-v7.3');
    
    % Create SUCCESS indicator
    fid = fopen(fullfile(runopt.path.me_mna_dae_ic, 'SUCCESS'),'w');
    fclose(fid);
    
else
    error('this has not yet been implemented')
end

Y0_consistent = y0;
Yp0_consistent = yp0;

end

