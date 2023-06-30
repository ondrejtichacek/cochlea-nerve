function [Y0_consistent, Yp0_consistent, mechopt] = InitialConditions( ...
    tspan, ODEFUN_ORIG, dep_fun, size_info, solveropt, mechopt, mnaopt, runopt, opt, ...
    y0, Y0, Sources, Q, M, Mred, B0, jac_pattern, ...
    system_const_params, jacob_const_params, common_const_params, ...
    PDEDiscretize_extra_args, args)
%INITIALCONDITIONS
arguments
    tspan
    ODEFUN_ORIG
    dep_fun
    size_info
    solveropt
    mechopt
    mnaopt
    runopt
    opt
    y0
    Y0
    Sources
    Q
    M
    Mred
    B0
    jac_pattern
    system_const_params
    jacob_const_params
    common_const_params
    PDEDiscretize_extra_args
    
    args.max_tries (1,1) double = 10
end

% for debug only
% y0 = Y0;

% extract original OdeFun_t, put it back in the end
OdeFun_t = mechopt.OdeFun_t;
mechopt.OdeFun_t = @(t) 0;

mechopt.vohc_0 = NaN;

system_const_params_IC = system_const_params;
jacob_const_params_IC = jacob_const_params;
jacob_const_params_IC_0 = jacob_const_params;
%%

jacob_const_params_IC_0{1} = B0;

%% currently not used

% ODEFUN = ODEFUN_ORIG;

%% increase BM damping
if false
    error('Also edit jacob const params!')

    PDEDiscretize_extra_args = [PDEDiscretize_extra_args, {'damping_factor', 1000}];
            
    [~, ~, ~, ~, ~, ~, ~, ~, ~, damp] = ...
            mech.v2.PDEDiscretize(mechopt, PDEDiscretize_extra_args{:});
    
    m = size_info.mech.start;
    n = size_info.mech.end;
    
    N = mechopt.Numstacks;
    
    A0 = system_const_params_IC{1};
    A0_mech = A0(m:n, m:n);
    
    % J22
    p = N+1;
    q = 2*N;
    
    A0_mech(p:q,p:q) = A0_mech(p:q,p:q) - (-damp);
    
    % J42
    p = N+1;
    q = 2*N;
    r = 3*N+1;
    s = 4*N;
    
    TMa = ones(N,1);
    A0_mech(r:s,p:q) = A0_mech(r:s,p:q) - (diag(TMa)*damp);
    
    A0(m:n, m:n) = A0_mech;
    
    system_const_params_IC{1} = A0;
end
%%

switch mechopt.integration
    case 'standalone'
        % Assembling f(Y,t), the RHS of the MNA equation
        ODEFUN = @(t, y, Z, t_mech, BM_displ, OHC_cilia) oc_mna.System(t, y, Z, ...
                        system_const_params_IC{:}, common_const_params{:}, ...
                        t_mech, BM_displ, OHC_cilia);

        % Caclulating the Jacobian term J = df/dy
        JACOBFUN = @(t, y, Z, t_mech, BM_displ, OHC_cilia) oc_mna.Jacob(t, y, Z, ...
                        jacob_const_params_IC{:}, common_const_params{:}, ...
                        t_mech, BM_displ, OHC_cilia);

        JACOBFUN_0 = @(t, y, Z, t_mech, BM_displ, OHC_cilia) oc_mna.Jacob(t, y, Z, ...
                        jacob_const_params_IC_0{:}, common_const_params{:}, ...
                        t_mech, BM_displ, OHC_cilia);

    case 'electric'
        
        % Assembling f(Y,t), the RHS of the MNA equation
        ODEFUN = @(t, y, Z) oc_mna.System(t, y, Z, ...
                system_const_params_IC{:}, common_const_params{:});

        % Caclulating the Jacobian term J = df/dy
        JACOBFUN = @(t, y, Z) oc_mna.Jacob(t, y, Z, ...
                jacob_const_params_IC{:}, common_const_params{:});
        
        JACOBFUN_0 = @(t, y, Z) oc_mna.Jacob(t, y, Z, ...
                jacob_const_params_IC_0{:}, common_const_params{:});
        
end

solveropt_copy = odeset(solveropt, ...
    'RelTol', 1e-3);

if ~isempty(solveropt_copy.Jacobian)
    solveropt_copy = odeset(solveropt_copy, ...
        'Jacobian', JACOBFUN_0);
end

% solveropt_copy = odeset( solveropt_copy, ...
%         'JPattern', jac_pattern);

%%

if true % currently not used
% This workaround is to fix an error of inconsistent initial condition:
% The DAE (MNA) itself does not have consistent Y0, while the ODE (channels)
% is consistent. This gives an error because the Y0_ODE is 'too good', and
% daeic3 fails to make it better. Therefore, we 'break' it below:

    if ~isempty(mnaopt.channels)
        a = size_info.channels.start;
        b = size_info.channels.end;

        Y0(a:b) = Y0(a:b)*(0.99);
        % Y0((mm-nn+1):end) = YYwrong;
    end
end

% [runopt, ~] = findResult('oc_mna_dae_ic', ...
%         [], [], [], mechopt, mnaopt, [], [], runopt, opt);
% 
% if runopt.found.oc_mna_ic == false
%     [runopt_tmp, ~] = findResult('oc_mna_dae_ic', ...
%         [], [], [], mechopt, mnaopt, [], [], copy(runopt), opt, ...
%         'result_base_dir', runopt.read_only_cache);
%     if runopt_tmp.found.mna_ic == true
%         runopt.found.oc_mna_ic = runopt_tmp.found.mna_ic;
%         runopt.path.oc_mna_ic = runopt_tmp.path.mna_ic;
%     end
% end

extra_dirs = {runopt.read_only_cache};

[runopt, ~] = findResultExtended(extra_dirs, 'oc_mna_dae_ic', ...
        [], [], [], mechopt, mnaopt, [], [], runopt, opt);

cache_dae_file = fullfile(runopt.path.oc_mna_dae_ic, 'dae_ic.mat');
if ~exist(runopt.path.oc_mna_dae_ic, 'dir')
    mkdir(runopt.path.oc_mna_dae_ic)
end

% Compute consistent inital conditions of the DAE for the solver
% if not continuing an existing simulation
if ~runopt.recalculate.oc_mna_dae_ic ...
        && isempty(y0) ...
        && runopt.found.oc_mna_dae_ic ...
        && exist(cache_dae_file, 'file')
    
    dbprintf('Loading consistent initial conditions from a file:\n -- %s.\n', cache_dae_file);
    tmp = load(cache_dae_file);
    y0 = tmp.y0;
    yp0 = tmp.yp0;
    
elseif runopt.no_create_only_load_oc_mna_dae_ic
    
    error('We require initial conditions to exist.')
    
elseif isempty(y0)

    OutputFcnList = {};

    if runopt.draw_gui == true
        
        audiotime = Time([0,10], 'ms');
        audio = zeros(size(audiotime));
        acc = zeros(size(audiotime));
    
        sig = [audio; acc];
    
        mechguiopt = mech.gui_v2.init( ...
            mnaopt.xgrid, audiotime, sig, mnaopt);
        
        mechguiopt.size_info = size_info;
        mechguiopt.abs_tol = Inf;
    
        mechguiopt.draw_period = Time(0, 'ms');
        
        OutputFcnList{end+1} = @(t,y,flag) mech.gui_v2.update(t, y, flag, ...
            mechguiopt, mechopt, mnaopt);
    end

   
    dbprintf('Computing consistent initial conditions.\n');
    timerInitialConditions = tic;
    
    try_again = true;
    try_num = 0;
    ic_system_method = 'M_fast';
    % ic_system_method = 'M_red_fast';
    % ic_system_method = 'M';
    % ic_system_method = 'M_full';
    % ic_system_method = 'no_M';
    
    % figure
    % hold on
    
    while try_again
        try_num = try_num + 1;
        
        dbprintf(' -- Using ''%s'' modification of the system to compute initial conditions\n', ic_system_method)
        
        switch ic_system_method
            case 'M_fast'
                Mfast = M;
                mm = size_info.circuit.start;
                nn = size_info.circuit.end;
                Mfast(mm:nn,mm:nn) = Mfast(mm:nn,mm:nn) * 1e-4;
                solveropt_init = odeset(solveropt_copy, 'Mass', Mfast); % 1-e4 factor to fasten the dynamics
            case 'M_red_fast'
                solveropt_init = odeset(solveropt_copy, 'Mass', Mred * 1e-4); % 1-e4 factor to fasten the dynamics
            case 'M'
                solveropt_init = odeset(solveropt_copy, 'Mass', M);
            case 'M_full'
                solveropt_init = solveropt_copy;
            case 'no_M'
                % Remove capacitance matrix because we don't need it to calculate consistent inital conditions.
                % solveropt_init = odeset(solveropt, 'Mass', []);
                rethrow(E)
            otherwise
                rethrow(E)
        end

        if try_num == 1
            tspan = Time([0, 10], 'us'); % super short tspan allows IC to pass daeic3
        elseif try_num < 10
            tspan = Time([0, 10], 'ms');
        else
            tspan = Time([0, 10], 'ms');
        end

        dbprintf(' -- tspan: [%s, %s]\n', tspan(1), tspan(2))

        tsp = tspan.s;
        
        OutputFcnList_init = OutputFcnList;

        if runopt.verbose >= 3
        % Display solver statistics
        % Waitbar function as OutputFcn, to be used after every integration timestep
        
            ws = WaitbarStorage();
            ws.num_div = mnaopt.NumDiv;
        
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

        switch mechopt.integration
            case 'standalone'
                % tmp = oc_mna.JacobLoadBM(0, fMechResult, unit_t, x, mnaopt.xgrid);
                % BM_displ = tmp{2};
                % OHC_cilia = tmp{3};
                BM_displ = zeros(1, mnaopt.Numstacks);
                OHC_cilia = zeros(1, mnaopt.Numstacks);
                extra_vars = {tsp, repmat(BM_displ, N, 1), repmat(OHC_cilia, N, 1)};
            case 'electric'
                extra_vars = {};
        end
        
        Y0_tmp = Y0;
            
        if try_num > 1
            % Extract IHC & OHC receptor potentials from solution
            fn = {'vihc', 'vohc'};
            for ii = 1:numel(fn)
                dep = fn{ii};
                vars.(dep) = dep_fun.(dep)(Y0(Q));                
            end

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

        % plot(Y0)
        % drawnow
        
        timerInitialConditions_iter = tic;

        % try
        sol = mtode.ode15s(ODEFUN, tsp, Y0_tmp, ...
                mtode.odeset(solveropt_init,  ...
                    'Stats', 'off'), ...
                Sources, extra_vars{:});
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
            
            if false
                if strcmp(fn{k}, 'mech')
                    figure
                    hold on
                    plot(Y0(m:n))
                    plot(y0(m:n))
                    drawnow
                end
            end
        end
        Y0 = y0;

        dbprintf(' -- iteration %d runtime %s.\n', try_num, disp_toc(toc(timerInitialConditions_iter)));
        
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
    
    max_expected_values = struct( ...
        'mech_BMx', 1e-9, ...
        'mech_BMv', 1e-6, ...
        'mech_TMx', 1e-9, ...
        'mech_TMv', 1e-6);
    
    fn = fieldnames(max_expected_values);
    for k = 1:numel(fn)
    
        m = size_info.(fn{k}).start;
        n = size_info.(fn{k}).end;
    
        if max(abs(Y0(m:n))) > max_expected_values.(fn{k})
            error('Initial conditions did not converge! Maximal expected value of %g for %s exceeded. Value = %g', ...
                max_expected_values.(fn{k}), fn{k}, max(abs(Y0(m:n))))
        end
    end
    
    %% Finalizing
    
    dbprintf('... runtime %s.\n', disp_toc(toc(timerInitialConditions)));
    % figure; surf(getVoltage(yy, 'IHC', mnaopt))
    
    fprintf('Saving initial guess for MNA DAE to %s\n', cache_dae_file)
    save(cache_dae_file, 'y0', 'yp0', '-v7.3');
    
    % Create SUCCESS indicator
    fid = fopen(fullfile(runopt.path.oc_mna_dae_ic, 'SUCCESS'),'w');
    fclose(fid);
    
else
    error('this has not yet been implemented')
end

Y0_consistent = y0;
Yp0_consistent = yp0;

m = size_info.mech.start;
n = size_info.mech.end;
Y0_consistent(m:n) = 0;

% Extract IHC & OHC receptor potentials from solution
fn = {'vihc', 'vohc'};
for ii = 1:numel(fn)
    dep = fn{ii};
    vars.(dep) = dep_fun.(dep)(Y0_consistent(Q));
end
mechopt.vohc_0 = vars.vohc;

mechopt.OdeFun_t = OdeFun_t;

end

% function calculate_error_norm()
% end