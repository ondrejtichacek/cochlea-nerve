function [ runopt ] = mainMNA( stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, y0 )
arguments
    stimulus (1,1) {mustBeA(stimulus, ["Signal", "CompoundSignal"])}
    topt (1,1) timeOpt
    midopt (1,1) midOpt
    mechopt (1,1) mechOpt
    mnaopt (1,1) mnaOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt
    y0 double = []
end
%MNA Modified Nodal Analysis of Cochlear RC network with mechanical input
%
% Solves the DAE:
%   M*dY/dt + A*Y - Z = 0.
%
% sometimes written as
%   M*dY/dt = -A*Y + Z =: f(Y,t)
%
% the matrix -A is often called the Jacobian, because
%   J := partial (f(Y,t)) / partial x = -A(t)
%
% However, when A = A(y,t) the product rule must be applied
%   J = -A(t,y) - dA(t,y)/dy * y
%
% In the above equations:
% * M : constant, singular mass matrix -- matrix of capacitances,
%       sometimes denoted by C
% * Y : vector of voltages at the nodes
% * A : block matrix of the form
%       |G B|
%       |E D|
%       * G : time dependent conductances, G(t) = G0 + delta_G(t)
%       * B,E : connectivity to voltage sources
%       * D : zero if only independent voltage sources are present
%
% * Z : vector of voltage and current sources
%       sometimes denoted by Sources
% 
% IMPLEMENTATION NOTES:
% * B,D,E,G and M are calculated outside the ode solver. Only delta_G(t)
%   needs to updated and added to block G0 of A0 to create A
% * OHC unstable with X1
%
%   See also ONEMNAJOB
%
% Copyright 2013-2021 OT
%

if isempty(opt.scratchdir)
    use_scratch = false;
else
    use_scratch = true;
end

unit_t = mnaopt.simulation_units.time;
unit_f = Time.inverse_unit(unit_t);

%% Create output time samples

% mnaopt.samplingFrequency = 64 * stimulus.frequency;

if mnaopt.evalAtTimePoints
    % create equidistant output samples
    numSteps = ceil((topt.compute.tf - topt.compute.t0) * mnaopt.samplingFrequency);
    numSamples = numSteps + 1;

    tspan = linspace(topt.compute.t0.(unit_t), topt.compute.tf.(unit_t), numSamples);
    tspan = tspan(:);
else
    tspan = topt.compute.tspan;
end

%% Mechanical model
% Calculate/load the mechanical results to be used by the electrical model.
% We solve the ode
%   M * dy/dt = J * y + g(t) + h(y)
% where g does not depend on y and and h is nonlinear function of y
% or the ode
%   M * dy/dt = J * y + g(t) + h(v)
% where h is a nonlinear function of voltage.
%

switch mechopt.integration
    case 'standalone'
        %  This version does not include feedback to the mechanical model
        
        disp_title( 'mechanical model' );
        [runopt, fMechResult, x] = MechResult( ...
                stimulus, copy(topt), midopt, mechopt, ...
                runopt, opt, memopt, paropt );

        stimulus_audio = stimulus;

        stimulus_stapes = MiddleEar(stimulus_audio, mechopt.middle_ear_params.identifier, ...
            topt.compute, topt, midopt, mechopt, runopt, opt, memopt, paropt);

        PDEDiscretize_extra_args = {};

    case 'electric'
        % Mechanical model incorporated to the electrical

        stimulus_audio = stimulus;

        stimulus_stapes = MiddleEar( ...
            stimulus_audio, mechopt.middle_ear_params.identifier, ...
            mechopt.middle_ear_params, mechopt.outer_ear_params, ...
            topt.compute, topt, midopt, runopt, opt, memopt, paropt);

        use_mna_xgrid = true;
        % use_mna_xgrid = false;
        if use_mna_xgrid == true
            PDEDiscretize_extra_args = {'xgrid', mnaopt.xgrid};
            assert(mechopt.Numstacks == mnaopt.Numstacks);
        else
            PDEDiscretize_extra_args = {};
        end
        
        mech_Numstacks = mechopt.Numstacks;

        time_unit = 's';
        stimulus_fcn = stimulus_stapes.eval_fcn(time_unit);

        [mechopt.xgrid, BMx, BMv, BMy, in, GG, TMy, TMw, TMa, damp, BMm] = ...
                mech.v2.PDEDiscretize(mechopt, ...
                    ... 'plotflag', true, ...
                    PDEDiscretize_extra_args{:});

        I = eye(mechopt.Numstacks);

        mechopt.mass = sparse(blkdiag(I, GG, I, GG));
%         mechopt.mass = sparse(blkdiag(I, GG, I, BMm));

        % this is actually a constant part of the jacobian
        % independent of the solution y
        mechopt.jacobian = sparse(mech.v2.Jacobian( ...
                mechopt.Numstacks, GG, BMx, BMv, BMy, TMy, TMw, TMa, ...
                mechopt.approximation, mechopt.amplifier));

        YCutDer = mechopt.YCutFunDer;

        prestin_voltage_sensitivity = mechopt.prestin_voltage_sensitivity;

        num_mech_variables = 4*mechopt.Numstacks;

        size_mech = num_mech_variables;

        Y0mech = zeros(num_mech_variables, 1);
        % Y0mech = 1e-6*ones(num_variables, 1);

        YCut = mechopt.YCutFun;

        mechopt.OdeFun_t = @(t) mech.v2.drivingForce( ...
            t, mech_Numstacks, in, TMa, [], stimulus_fcn, []);

        mechopt.OdeFun_y = @(y) mech.v2.amplifierNonlinearity( ...
            [], y, mech_Numstacks, [], TMa, BMy, [], YCut);

        mechopt.jacobian_nonlin_part = @(y) mech.v2.JacobianNonlinPart( ...
            y, mech_Numstacks, BMy, TMa, YCutDer);

        if mechopt.amplifier == "electric"
            if mechopt.approximation == "linear"
                mechopt.OdeFun_v = @(V) mech.v2.amplifier_lin_VOHC( ...
                    [], V, mech_Numstacks, [], TMa, BMy, [], [], prestin_voltage_sensitivity);

                mechopt.jacobian_lin_part_vohc = @() mech.v2.Jacobian_lin_VOHC( ...
                    BMy, TMa, prestin_voltage_sensitivity);
                
            elseif mechopt.approximation == "nonlinear"
                mechopt.OdeFun_v = @(V) mech.v2.amplifierNonlinearity_VOHC( ...
                    [], V, mech_Numstacks, [], TMa, BMy, [], YCut, prestin_voltage_sensitivity);

                mechopt.jacobian_nonlin_part_vohc = @(V) mech.v2.JacobianNonlinPart_VOHC( ...
                    V, mech_Numstacks, BMy, TMa, YCutDer, prestin_voltage_sensitivity);

            else
                error('Unknown option %s', mechopt.approximation)
            end
        end

        debug_try = false;
        % debug_try = true;

        if debug_try

            maxstep = 1/Frequency(44.1, 'kHz');
            initialstep = 1/Frequency(200, 'kHz');

            solveropt = odeset( ...
                'RelTol', 1e-6, ...
                'AbsTol', mech_err_tol(mechopt, stimulus), ...
                'MaxStep', maxstep.(unit_t), ...
                'InitialStep', initialstep.(unit_t));

            J = mechopt.jacobian;
            g = mechopt.OdeFun_t;
            h = mechopt.OdeFun_y;

            MECHODE = @(t,y) J * y + g(t) + h(y);

            % N = 30;
            % tsp = tspan(1:N);
            tsp = tspan;

            solveropt = odeset( solveropt, ...
                'BDF', 'on', ...
                'MaxOrder', 2, ...
                'Mass', mechopt.mass, ...
                'Jacobian', @(t,y) mechopt.jacobian + mechopt.jacobian_nonlin_part(y), ...
                'Stats','off' ...
                );

            sig = [stimulus.audio(:)];

            mechguiopt = mech.gui.init(mechopt.xgrid, mechopt.Numstacks, stimulus.audiotime, sig, NaN, [], {'BM', 'TM'});
            OutputFcn = @(t, y, flag) mech.gui_v1.update(t, y, flag, mechopt.Numstacks, mechopt.NMkonst_BM, mechguiopt);

            solver = @ode15s;
            % solver = @ode45;
            
            sol = solver(MECHODE, tsp, Y0mech, ...
                    odeset(solveropt, 'OutputFcn', OutputFcn, 'Stats', 'off'));

            figure, mesh(mechopt.NMkonst_BM*sol.y(1:mechopt.Numstacks,:))
            %figure, mesh(mechopt.NMkonst_BM*sol.y(mechopt.Numstacks+1:2*mechopt.Numstacks,:))
            %figure, mesh(mechopt.NMkonst_OHC_cilia*sol.y(2*mechopt.Numstacks+1:3*mechopt.Numstacks,:))
            %figure, mesh(mechopt.NMkonst_OHC_cilia*sol.y(3*mechopt.Numstacks+1:4*mechopt.Numstacks,:))

        end
        
    otherwise
        error('Unsupported value of mechopt.integration = %s', mechopt.integration)
end

%% Electrical model
disp_title( 'electrical model' );

%% MNA initialization
% M*dY/dt = J*Y + Sources
% Calculate the constant M, A0 (A without delta_G(t)), Sources and Y0, the
% initial voltages.

dbprintf('Assembling circuit\n');
timerCircuitAssembly = tic;

extra_dirs = {runopt.read_only_cache};

[runopt, ~] = findResultExtended(extra_dirs, 'oc_mna_circuit_ic', ...
        [], [], [], mechopt, mnaopt, [], [], runopt, opt);

cache_circuit_file = fullfile(runopt.path.oc_mna_circuit_ic, 'circuit_ic.mat');
if ~exist(runopt.path.oc_mna_circuit_ic, 'dir')
    mkdir(runopt.path.oc_mna_circuit_ic)
end

if ~runopt.recalculate.oc_mna_circuit_ic ...
        && runopt.found.oc_mna_circuit_ic ...
        && exist(cache_circuit_file, 'file')
    dbprintf('Loading circuit voltage initial guess from a file:\n -- %s.\n', cache_circuit_file);
    IC = load(cache_circuit_file);
elseif runopt.no_create_only_load_oc_mna_circuit_ic
    
    error('We require initial conditions to exist.')
    
else
    IC = struct();
end

[M, M0, ~, A0, Sources, Y0, p, q, vihc, vohc] = circuit_equations( ...
        mnaopt, 'initial_guess', IC);

if ~isempty(fieldnames(IC))
    dbprintf('Norm of initial guess changed by:\n')
    dbprintf('    vihc: %g\n', norm(IC.vihc - vihc) / norm(IC.vihc));
    dbprintf('    vohc: %g\n', norm(IC.vohc - vohc) / norm(IC.vohc));
else
    dbprintf('Saving initial guess for circuit_equations to %s\n', cache_circuit_file)
    save(cache_circuit_file, 'vihc', 'vohc', '-v7.3');
    % Create SUCCESS indicator
    fid = fopen(fullfile(runopt.path.oc_mna_circuit_ic, 'SUCCESS'),'w');
    fclose(fid);
end
    
dbprintf('... runtime %s.\n', disp_toc(toc(timerCircuitAssembly)));

size_MNA = numel(Y0);

spec_size = mnaopt.circuit.num_v + mnaopt.circuit.num_node;

size_info.circuit = struct( ...
    'spec_size', spec_size, ...
    'size', size_MNA, ...
    'start', 1, ...
    'end', size_MNA);

size_info.blocks = "circuit";

%% Solution of DAE

disp_title( 'solving DAE', '-' );

if use_scratch
    DIR = strrep(runopt.path.oc_mna, opt.resDir, opt.scratchdir);
    if exist(DIR, 'dir')
        srmdir(DIR, 's');
    end
    mkdir(DIR);
else
    DIR = runopt.path.oc_mna;
end

res.y = matfile(fullfile(DIR,'Y.mat'),'Writable',true);
res.t = matfile(fullfile(DIR,'t.mat'),'Writable',true);

checkpoints.y = matfile(fullfile(DIR,'cp_Y.mat'),'Writable',true);
checkpoints.t = matfile(fullfile(DIR,'cp_t.mat'),'Writable',true);

%% Set up problem functions

xx = mnaopt.xgrid(:);

% Jacobian initialization:
% precalculate functions for constructing A from A0 (by adding delta_G(t))
[g, dg, Ig, c, Ic, dep_fun, dep_sel, dep_sign] ...
        = oc_mna.JacobPrecalc(mnaopt.circuit, mnaopt.NumstacksScaling, ...
                mnaopt.xgrid, mnaopt);

% inverse permutations
P(p) = 1:numel(p);
Q(q) = 1:numel(q);

% Extract IHC & OHC receptor potentials from solution
fn = {'vihc', 'vohc'};
for ii = 1:numel(fn)
    dep = fn{ii};
    vars.(dep) = dep_fun.(dep)(Y0(Q));
end

% Link to steady-state solution
vars.vihc_ss = vars.vihc;
vars.vohc_ss = vars.vohc;

%% Extension -- basolateral IHC channels open probability as 2nd order ODE

channels = mnaopt.channels;

if ~isempty(channels)
    V = struct('IHC', vars.vihc, 'OHC', vars.vohc, 'IHC_ss', vars.vihc_ss, 'OHC_ss', vars.vohc_ss);
    [MM, AA, zz] = ChannelspOpen(V, channels, mnaopt.Numstacks);

    size_channels = numel(zz);

    [a,b] = size(M);
    [aa,bb] = size(MM);

    Mred = [         M, sparse(a,bb);
         sparse(aa,b), sparse(aa,bb)];

    M =  [         M, sparse(a,bb);
         sparse(aa,b),          MM];

    M0 = [        M0, zeros(a,bb);
         zeros(aa,b), spalloc(aa,bb,nnz(MM))];

    A0 = [        A0, spalloc(a,bb, size_channels);
         zeros(aa,b), spalloc(aa,bb,nnz(AA))];

    YY = AA\zz;
    % YYwrong = AA\-zz;

    dbprintf('Initial condition found with residue norm = %f\n', norm(AA*YY-zz))

    Y0 = [Y0; YY];

    assert(aa == mnaopt.Numstacks*sum([channels.order]))

    % update ordering
    p = [p, numel(p) + (1:aa)];
    q = [q, numel(q) + (1:aa)];

    s = 0;
    for block = size_info.blocks
        s = s + size_info.(block).size;
    end
    
    size_info.channels = struct( ...
        'size', size_channels, ...
        'start', s + 1, ...
        'end', s + size_channels);
    
    size_info.blocks = [size_info.blocks, "channels"];

    dep_names{numel(channels)} = [];
    dep_names_alt{numel(channels)} = [];

    channel_names = {channels.name};
    var_names = {channels.var_name};
    for k = 1:numel(var_names)
        cn = channel_names{k};
        v = var_names{k};
        
        vv = ['channel_', cn];
        vn = ['popen_', cn];

        dep_names{k} = vn;
        dep_names_alt{k} = vv;

        n = 2 * mechopt.Numstacks;
            
        size_info.(vv) = struct( ...
            'size', n, ...
            'start', size_info.channels.start + (k-1)*n , ...
            'end', size_info.channels.start + (k)*n - 1);

        o = 1; % solution O ... open propability is in the first equation
        % o = 2; if the solution was always in the second equation

        dep_sel.(v) = (size_info.(vv).start + (o - 1)) : 2 : size_info.(vv).end;
    end
end

%% Init mass pattern

% Boolean pattern matrix: 1 when M(y,t) depends on y
mass_pattern = oc_mna.MassPattern(Ic, channels, mnaopt.Numstacks, spec_size);


%% Extend system by F = FOHC ~ high pass filtered VOHC
%
% RC dF/dt + F = RC dV/dt
%
%    F ... OHC force
%    V ... VOHC = V^e-V^i

if true%false

    s = 0;
    for block = size_info.blocks
        s = s + size_info.(block).size;
    end
    
    dVdt = sparse(1:mnaopt.Numstacks, dep_sel.vohc{1}, dep_sign.vohc(1), mnaopt.Numstacks, s) ...
        + sparse(1:mnaopt.Numstacks, dep_sel.vohc{2}, dep_sign.vohc(2), mnaopt.Numstacks, s);
    
    dFdt = speye(mnaopt.Numstacks, mnaopt.Numstacks);
    
    tau = 20;
    RC = 1/tau;
    
    AA = speye(mnaopt.Numstacks, mnaopt.Numstacks);
    
    YY0 = zeros(mnaopt.Numstacks, 1);
    
    [a,b] = size(M);
    [aa,bb] = size(AA);
    
    M =  [         M, sparse(a,bb);
            -RC*dVdt, RC*dFdt];
    %           -dVdt,     dFdt];
    
    Mred = M;
    
    M0 =  [       M0, sparse(a,bb);
            -RC*dVdt, RC*dFdt];
    %            -dVdt,     dFdt];
    
    if ~isempty(mass_pattern)
        mass_pattern = sparse(blkdiag(mass_pattern, spalloc(mnaopt.Numstacks,mnaopt.Numstacks,0)));
    end
    
    A0 = [        A0, sparse(a,bb);
         sparse(aa,b), AA];
    
    Y0 = [Y0; YY0];
    
    % update ordering
    p = [p, numel(p) + (1:mnaopt.Numstacks)];
    q = [q, numel(q) + (1:mnaopt.Numstacks)];
    
    size_dvohc = mnaopt.Numstacks;
    
    size_info.dvohc = struct( ...
        'size', size_dvohc, ...
        'start', s + 1, ...
        'end', s + size_dvohc);
    
    size_info.blocks = [size_info.blocks, "dvohc"];
end

mechopt.vihc_ss = vars.vihc;
mechopt.vohc_ss = vars.vohc;

%% Extend system by mechanical PDE

switch mechopt.integration
    case 'standalone'
    case 'electric'
        mechopt.vohc_0 = vars.vohc;

        MM = mechopt.mass;

        [a,b] = size(M);
        [aa,bb] = size(MM);

        Mred = [         M, sparse(a,bb);
             sparse(aa,b), sparse(aa,bb)];

        M =  [         M, sparse(a,bb);
             sparse(aa,b),          MM];

        M0 = [        M0, zeros(a,bb);
             zeros(aa,b),           MM];

        if ~isempty(mass_pattern)
            mass_pattern = sparse(blkdiag(mass_pattern, spalloc(num_mech_variables,num_mech_variables,0)));
        end

        % constant part
        AA = - mechopt.jacobian;
        % AA = mechopt.jacobian;

        nz = 2 * 4 * mechopt.Numstacks;
        
        A0 = [        A0, spalloc(a,bb, nz);
             spalloc(aa,b, nz), AA];

        Y0 = [Y0; Y0mech];

        % update ordering
        p = [p, numel(p) + (1:num_mech_variables)];
        q = [q, numel(q) + (1:num_mech_variables)];

        s = 0;
        for block = size_info.blocks
            s = s + size_info.(block).size;
        end

        size_info.mech = struct( ...
            'size', size_mech, ...
            'start', s + 1, ...
            'end', s + size_mech);

        size_info.blocks = [size_info.blocks, "mech"];

        var_names = {'BMx', 'BMv', 'TMx', 'TMv'};
        alt_var_names = {'BM_displ', 'BM_velocity', 'OHC_cilia', 'OHC_cilia_velocity'};
        for k = 1:numel(var_names)
            v = var_names{k};
            vv = ['mech_', v];
            size_info.(vv) = struct( ...
                'size', mechopt.Numstacks, ...
                'start', s + 1 + (k-1)*mechopt.Numstacks, ...
                'end', s + (k)*mechopt.Numstacks);
            
            w = alt_var_names{k};
            dep_sel.(w) = size_info.(vv).start : size_info.(vv).end;
        end
        
    otherwise
        error('Unsupported value of mechopt.integration = %s', mechopt.integration)
end

%% Extend system by IHC_stereocilia S = high pass filtered BMx
%
% RC dS/dt + S = RC dX/dt
%
%    S ... IHC_stereocilia ... high pass filtered BMx
%    X ... BMx

if false

    s = 0;
    for block = size_info.blocks
        s = s + size_info.(block).size;
    end
    
    ss = size_info.mech_BMx.start : size_info.mech_BMx.end;

    dVdt = sparse(1:mechopt.Numstacks, ss, 1, mechopt.Numstacks, s);
    
    dFdt = speye(mechopt.Numstacks, mechopt.Numstacks);
    
    tau = 20;
    RC = 1/tau;
    
    AA = speye(mechopt.Numstacks, mechopt.Numstacks);
    
    YY0 = zeros(mechopt.Numstacks, 1);
    
    [a,b] = size(M);
    [aa,bb] = size(AA);
    
    M =  [         M, sparse(a,bb);
            -RC*dVdt, RC*dFdt];
    %           -dVdt,     dFdt];
    
    Mred = M;
    
    M0 =  [       M0, sparse(a,bb);
            -RC*dVdt, RC*dFdt];
    %            -dVdt,     dFdt];
    
    if ~isempty(mass_pattern)
        mass_pattern = sparse(blkdiag(mass_pattern, spalloc(mechopt.Numstacks,mechopt.Numstacks,0)));
    end
    
    A0 = [        A0, sparse(a,bb);
         sparse(aa,b), AA];
    
    Y0 = [Y0; YY0];
    
    % update ordering
    p = [p, numel(p) + (1:mechopt.Numstacks)];
    q = [q, numel(q) + (1:mechopt.Numstacks)];
    
    size_IHC_stereocilia = mechopt.Numstacks;
    
    size_info.IHC_stereocilia = struct( ...
        'size', size_IHC_stereocilia, ...
        'start', s + 1, ...
        'end', s + size_IHC_stereocilia);
    
    size_info.blocks = [size_info.blocks, "IHC_stereocilia"];
end

%% Extend system by regularization term

if true%false

    s = 0;
    for block = size_info.blocks
        s = s + size_info.(block).size;
    end
    
    YY0 = zeros(1, 1);
    
    AA = 0;

    [a,b] = size(M);
    [aa,bb] = size(AA);

    M =  [           M, sparse(a,bb);
           sparse(aa,b),  0];
    
    Mred = M;
    
    M0 =  [ M0, sparse(a,bb);
            sparse(aa,b),  0];
    
    if ~isempty(mass_pattern)
        mass_pattern = sparse(blkdiag(mass_pattern, spalloc(1,1,0)));
    end
    
    A0 = [        A0, sparse(a,bb);
         sparse(aa,b), AA];
    
    Y0 = [Y0; YY0];
    
    % update ordering
    p = [p, numel(p) + (1:1)];
    q = [q, numel(q) + (1:1)];
    
    size_reg = 1;
    
    size_info.reg = struct( ...
        'size', size_reg, ...
        'start', s + 1, ...
        'end', s + size_reg);
    
    size_info.blocks = [size_info.blocks, "reg"];
end

%%

% since the system is solved in the form
%    C*dy/dy = -A(t,y) + z(t,y)
% let us denote
%    -A(t,y) = B(t,y)

B0 = spalloc(size(A0,1), size(A0,2), nzmax(A0));
B0(1:end, 1:end) = -A0;

% similar thing for the constant part of the jacobian

jac_0 = oc_mna.Jacob_const(mnaopt.Numstacks, size_info, dep_sel, dep_sign, p, q, mechopt);

nzmax_extra = numel(channels) * 2 * 2 * mnaopt.Numstacks;

J0 = spalloc(size(A0,1), size(A0,2), nzmax(A0) + nzmax_extra);
J0(1:end, 1:end) = -A0;
if ~isempty(jac_0)
    J0(1:end, 1:end) = J0(1:end, 1:end) + jac_0;
end

assert(nzmax(B0) == nzmax(A0))
assert(nzmax(J0) == nzmax(A0) + nzmax_extra)

system_const_params = {B0};
jacob_const_params = {J0};

if all(p == 1:numel(p)) && all(q == 1:numel(q))
    needs_ordering = false;
else
    needs_ordering = true;
end

% Constant parameters used to construct B from B0
common_const_params = {channels, mnaopt.Numstacks, size_info, mechopt.samplingFrequency.(unit_f), ...
                    g, dg, Ig, xx, ...
                    dep_fun, dep_sel, dep_sign, dep_names, dep_names_alt, ...
                    p, q, needs_ordering, mechopt};

% Constant parameters used to construct M from M0
mass_const_params = {M0, channels, mnaopt.Numstacks, spec_size, c, Ic, xx, ...
                   dep_fun, p, q, mechopt};
                
% Sparsity pattern of A, can be used by the solver to numerically calculate
% the nonzero elements of the Jacobian
% jac_pattern = oc_mna.JacobPattern(A0, channels, mnaopt.Numstacks, ...
%         g_BM_displ, I_BM_displ, g_OHC_cilia, I_OHC_cilia, ...
%         g_vihc, dg_vihc, I_vihc, g_vohc, dg_vohc, I_vohc, xx);

switch mechopt.integration
    case 'standalone'
        % Assembling f(Y,t), the RHS of the MNA equation
        ODEFUN = @(t, y, Z, t_mech, BM_displ, OHC_cilia) oc_mna.System(t, y, Z, ...
                        system_const_params{:}, common_const_params{:}, ...
                        t_mech, BM_displ, OHC_cilia);

        % Caclulating the Jacobian term J = df/dy
        JACOBFUN = @(t, y, Z, t_mech, BM_displ, OHC_cilia) oc_mna.Jacob(t, y, Z, ...
                        jacob_const_params{:}, common_const_params{:}, ...
                        t_mech, BM_displ, OHC_cilia);
                   
        if mnaopt.add_solveropt.use_jpattern == true
            jac_pattern = JACOBFUN(0,Y0,Sources,0,zeros(size(mnaopt.xgrid)),zeros(size(mnaopt.xgrid)));
            jac_pattern = sparse(jac_pattern ~=0);
        end

    case 'electric'
        
        % Assembling f(Y,t), the RHS of the MNA equation
        ODEFUN = @(t, y, Z) oc_mna.System(t, y, Z, ...
            system_const_params{:}, common_const_params{:});

        % Caclulating the Jacobian term J = df/dy
        JACOBFUN = @(t, y, Z) oc_mna.Jacob(t, y, Z, ...
            jacob_const_params{:}, common_const_params{:});
        
        % jac_pattern = JACOBFUN(0,Y0,Sources);
        if true || mnaopt.add_solveropt.use_jpattern == true
            jac_pattern = JACOBFUN(0,Y0,Sources);
            jac_pattern = sparse(jac_pattern ~=0);
        end
end
            
MASSFUN = @(t, y, Z, t_mech, BM_displ, OHC_cilia) oc_mna.CapacitanceMatrix( t, y, Z, ...
                mass_const_params{:});

%% Configure solver

if isa(stimulus, 'PureTone')
    % maxstep = 1/stimulus.frequency / 128;
    % maxstep = 1/stimulus.frequency / 64;
    maxstep = 1/max([stimulus.frequency]) / 32;
elseif isa(stimulus, 'Click') || isa(stimulus, 'Click2')
    maxstep = 1/mnaopt.samplingFrequency;
else
    % maxstep = 1/Frequency(44.1, 'kHz');
    maxstep = 1/mnaopt.samplingFrequency;
end

initialstep = min([maxstep, 1/mnaopt.samplingFrequency;]);

pos_ind = [ ...
    circuit_pos_ind(mnaopt), ...
    channel_pos_ind(mnaopt)];

if true%false
    pos_ind = [ pos_ind, ...
        dvohc_pos_ind(mnaopt)];
end

if strcmp(mechopt.integration, 'electric')
    pos_ind = [ pos_ind, ...
        mech_pos_ind(mechopt)];
end

pos_ind = {pos_ind(:)};

system_size = 0;
for ii = 1:numel(size_info.blocks)
    block = size_info.blocks(ii);
    system_size = system_size + size_info.(block).size;
end

system_size = system_size - size_info.reg.size;

% Loc scale
if isa(stimulus, 'PureTone') && numel(stimulus) == 1

    loc_scale = location_scale( ...
        stimulus.frequency.Hz, stimulus.amplitude, mechopt.gain, ...
        mechopt.xgrid, 0.06, 'plotflag', true);
    
    loc_scale_t_full = @(t) ones(system_size, 1);

else
    loc_scale = ones(1, mechopt.Numstacks);

    [S,f,t] = spectrogram(stimulus.audio, Time(3, 'ms') * stimulus.fs, [], [], stimulus.fs.Hz);
    [~, ind] = max(abs(S));
    % figure
    % plot(t, f(ind))

    t = t(:);
    f = f(:);

    FF = griddedInterpolant(t, f(ind), 'linear', 'nearest');

    % loc scale base to allow higher tollerance at the apex -- needed for
    % high frequency tones
    loc_scale_base = location_scale(20e3, 0, mechopt.gain, mechopt.xgrid, 0.09) - 0.1;
    loc_scale_base = 100 * loc_scale_base;

    loc_scale_t = @(t) loc_scale_base + location_scale(FF(t), stimulus.amplitude, mechopt.gain, mechopt.xgrid, 0.06*10);
        
    S = struct('type', '()', 'subs', {pos_ind});
    loc_scale_t_full = @(t) subsref(loc_scale_t(t)', S);
end

loc_scale_demo = ones(1, mechopt.Numstacks);

% Absolute tolerance
abs_tol = [ ...
    circuit_err_tol(mnaopt, stimulus, loc_scale), ...
    channel_err_tol(mnaopt, stimulus, loc_scale)];

if true%false
    abs_tol = [ abs_tol, ...
        dvohc_err_tol(mnaopt, stimulus, loc_scale)];
end

if strcmp(mechopt.integration, 'electric')
    abs_tol = [ abs_tol, ...
        mech_err_tol(mechopt, stimulus, loc_scale)];
end

%% Set tolerance according to envelope

tolerance_rel_envelope = true;

if tolerance_rel_envelope == true

    F = envelope_tolerance(stimulus, tspan, ...
        mnaopt.simulation_units.time, ...
        mnaopt.samplingFrequency, ...
        'amplifier', mechopt.gain > 0, ...
        'plotflag', true);

    G = envelope_tolerance(stimulus, tspan, ...
        mnaopt.simulation_units.time, ...
        mnaopt.samplingFrequency, ...
        'amplifier', mechopt.gain > 0, ...
        'saturate_spl', 70, ...
        'plotflag', true);
    
    abs_tol = abs_tol * 1;

    abs_tol = abs_tol / 10*2;

    abs_tol = abs_tol(:); % IMPORTANT

%     abs_tol = @(t) F(t) * abs_tol;

    abs_tol = @(t) [
        G(t) * abs_tol(size_info.circuit.start:size_info.circuit.end)
        % loc_scale_t(t) * F(t) * abs_tol(size_info.channels.start:size_info.channels.end)
        abs_tol(size_info.channels.start:size_info.channels.end)
        F(t) * abs_tol(size_info.dvohc.start:size_info.dvohc.end)
        F(t) * abs_tol(size_info.mech.start:size_info.mech.end)
        0.1 % reg
        ] .* [loc_scale_t_full(t); 1];
end


%%

% abs_tol = inf;
% abs_tol = eps;

% default solver options, might be overridden by mnaopt.solveropt
solveropt = odeset( ...
    'MaxStep', maxstep.(unit_t), ...
    'InitialStep', initialstep.(unit_t), ...
    'AbsTol', abs_tol, ...
    'RelTol', 1e-3, ...
    'BDF', 'on', ...
    'MaxOrder', 2, ...
    ...'Vectorized', 'on', ...
    'Stats','off' ...
    );


% Update Mass Matrix properties

switch mnaopt.capacitance_state_dependence
    case 'none'
        MASS_odeopts = {
            'Mass', M, ...                        % Capacitances
            ...'Mass', D, ...
            'MassSingular', 'yes', ...            % Singular matrix -> DAE
            'MStateDependence', 'none', ...       % M != M(y)
        };
    case 'weak'
        
        if isempty(mass_pattern)
            error('mnaopt.capacitance_state_dependence is set to %s but mass pattern is empty', mnaopt.capacitance_state_dependence)
        end
        
        MASS_odeopts = {
            'Mass', MASSFUN, ...                  % Capacitances
            'MassSingular', 'yes', ...            % Singular matrix -> DAE
            'MStateDependence', 'weak', ...       % M = M(y)
        };
    case 'strong'
        
        if isempty(mass_pattern)
            error('mnaopt.capacitance_state_dependence is set to %s but mass pattern is empty', mnaopt.capacitance_state_dependence)
        end
        
        MASS_odeopts = {
            'Mass', MASSFUN, ...                  % Capacitances
            'MassSingular', 'yes', ...            % Singular matrix -> DAE
            'MStateDependence', 'strong', ...     % M = M(y)
            'MvPattern', mass_pattern, ...
        };
    otherwise
        error('Unrecognised option %s', mnaopt.capacitance_state_dependence);
end

solveropt = odeset(solveropt, MASS_odeopts{:});

if mnaopt.zero_capacitance_test == true
    warning('Setting Capacitance matrix to zero -- use for testing only!')
    solveropt = odeset(solveropt, 'Mass', M * 1e-4);
end

% Update Jacobian properties

if mnaopt.add_solveropt.use_jacobian == true
    solveropt = odeset( solveropt, ...
        'Jacobian', JACOBFUN);
    dbprintf('Using analytic Jacobian.\n');
elseif mnaopt.add_solveropt.use_jpattern == true
    solveropt = odeset( solveropt, ...
        'JPattern', jac_pattern);
    dbprintf('Using numeric approximation of the Jacobian, employing the JPattern.\n');
else
    dbprintf('Using numeric approximation of the Jacobian.\n');
end

OutputFcnList = {};

if runopt.draw_gui == true
        
    audiotime = stimulus_stapes.audiotime;
    audio = stimulus_stapes.audio;
    acc = stimulus_stapes.eval(audiotime);

    if isempty(audiotime)
        audiotime = stimulus.audiotime;
        audio = stimulus.audio;
        acc = stimulus_stapes.eval(audiotime);
    end

    sig = [audio; acc];
    sig(2,:) = sig(2,:) / max(sig(2,:));

    mechguiopt = mech.gui_v2.init( ...
        mnaopt.xgrid, audiotime, sig, mnaopt, 'frequency', max([stimulus.frequency]));
    
    mechguiopt.size_info = size_info;
    mechguiopt.abs_tol = abs_tol;

    if isa(stimulus, 'PureTone')
        mechguiopt.draw_period = 1 / max([stimulus.frequency]) / 24;
    else
        mechguiopt.draw_period = Time(0.01, 'ms');
    end
    
    OutputFcnList{end+1} = @(t,y,flag) mech.gui_v2.update(t, y, flag, ...
        mechguiopt, mechopt, mnaopt);
end

if runopt.verbose >= 3
    % Display solver statistics
    solveropt = odeset(solveropt, 'Stats', 'on');

    ws = WaitbarStorage();
    ws.num_div = mnaopt.NumDiv;

    % Waitbar function as OutputFcn, to be used after every integration timestep
    if runopt.waitbarFunctionAvailable
        OutputFcnList{end+1} = @(t, y, flag) odewaitbar(t, y, flag, tspan(end), ws);
    else
        OutputFcnList{end+1} = @(t, y, flag) odewaittext(t, y, flag, tspan(end), ws);
    end
end

if numel(OutputFcnList) > 0
    OutputFcn = @(t,y,flag) OutputFcnChain(OutputFcnList, t, y, flag);
    
    % number of ~ must correspond with the number of extra arguments in ode15s
    OutputFcnWrapper = @(t, y, flag, ~, ~, ~, ~) OutputFcn(t, y, flag);
    solveropt = odeset(solveropt, 'OutputFcn', OutputFcnWrapper);
end

% overwrite default solver options by those specified in mnaopt
solveropt = odeset(solveropt, mnaopt.solveropt{:});

disp(solveropt)

%% Solve problem

wrapopt.append = false;
wrapopt.res = res;
wrapopt.verbose = runopt.verbose;
wrapopt.NumDiv = mnaopt.NumDiv;
switch mechopt.integration
    case 'standalone'
        wrapopt.InterDivFcn = @(tspan) oc_mna.JacobLoadBM(tspan, fMechResult, unit_t, x, mnaopt.xgrid);
    case 'electric'
end
wrapopt.memopt = memopt;
wrapopt.save_method = mnaopt.save_method; %'matlab_matfile';
wrapopt.save_checkpoints = true;
wrapopt.checkpoint_period = topt.compute.checkpoint_period.(unit_t);
wrapopt.checkpoints = checkpoints;
wrapopt.fallback_solver = mnaopt.fallback_solver;

% Get the length of saved checkpoints
if exist(checkpoints.t.Properties.Source, 'file')
    wrapopt.checkpoints_last_t_index = max(matfileWhosByName(checkpoints.t, 't', 'size')); % = memmory friendly length(checkpoints.t.t)
else
    wrapopt.checkpoints_last_t_index = 0;
end

own_initial_conditions = true;

if own_initial_conditions == true
    % Initial Conditions
    % note - mechopt is modified and used (it is a handle object)
    [Y0_consistent, Yp0_consistent, mechopt] = oc_mna.InitialConditions( ...
            tspan, ODEFUN, dep_fun, size_info, ...
            solveropt, mechopt, mnaopt, runopt, opt, ...
            y0, Y0, Sources, Q, M, Mred, B0, jac_pattern, ...
            system_const_params, jacob_const_params, ...
            common_const_params, PDEDiscretize_extra_args);
    
    % reg term fix.
    if numel(Y0_consistent) == numel(Y0) - 1
        Y0_consistent = [Y0_consistent; 0];
        Yp0_consistent = [Yp0_consistent; 0];
    end

    % Set initial slope to fix daeic3 inconsistent initial condition
    % M*yp0 = f
    solveropt = odeset(solveropt, 'InitialSlope', Yp0_consistent);
    
    solveropt = mtode.odeset(solveropt, 'SkipICCheck', true);
else
    Y0_consistent = Y0;
end
switch mnaopt.solver
    case {'ode15s', 'mtode.ode15s'}
        mnaopt.solver = @mtode.ode15s;
        % mnaopt.solver = @ode15s;
    case {'ode23t', 'mtode.ode23t'}
        mnaopt.solver = @mtode.ode23t;
        % mnaopt.solver = @ode23t;
    otherwise
        error('Unsupported solver %s', mnaopt.solver)
end

% Alternatively, 'break' the ODE part of the initial condition
% Y0_consistent((mm-nn+1):end) = YYwrong;

%% Freq. domain attempt

% [~, A, Z] = ODEFUN(0, Y0_consistent, Sources);
% M = MASSFUN(0, Y0_consistent, Sources, [], [], []);
% 
% Y = oc_mna.MNA_freq(Frequency(1000, 'Hz'), mnaopt.Numstacks, M, A, Z, TMa, in, size_info);
% 
% n = size_info.mech_BMx.start;
% m = size_info.mech_BMx.end;
% 
% YY = abs(Y(n:m));

%%

wrapopt.v2 = true;
if wrapopt.v2 == true
    wrapopt.maxt = 'auto';
    % wrapopt.maxt = 1000; % for testing purposes
    
    wrapopt.storage = odeWrapperStorage(tspan, Y0_consistent, wrapopt);

    if exist('OutputFcn', 'var')
        wrapopt.storage.output_functions{end+1} = OutputFcn;
    end

    OutputFcnWrapper = @(t, y, flag, ~, ~, ~, ~) wrapopt.storage.evalOutputFunctions(t, y, flag);
    solveropt = mtode.odeset(solveropt, 'OutputFcn', OutputFcnWrapper);
end


dbprintf('Solving the DAE.\n');
timerOdeWrapper = tic;

% Perform the actual calculation
% The results of the numerical integration are stored in wrapopt.storage.
% Example:
%   wrapopt.storage.ft.t : timesteps
%   wrapopt.storage.fy.y : voltages and currents for all the nodes in the MNA matrix
statistics = odeWrapper(wrapopt, mnaopt.solver, ...
        ODEFUN, tspan, Y0_consistent, solveropt, Sources);

dbprintf('\nTotal odeWrappper runtime %s.\n', disp_toc(toc(timerOdeWrapper)));
dbprintf('... ode runtime: %s.\n', disp_toc(sum([statistics.ode_time])));
dbprintf('... save time:   %s.\n', disp_toc(sum([statistics.save_time])));

% update runtime parameters
% runopt.update('datetime', clock, 'MNAtexec', texec);

%% Save parameters

res.y.ordering = struct('p', p, 'q', q);

res.y.size_info = size_info;

res.t.unit = unit_t;
checkpoints.t.unit = unit_t;

end_time = Time(res.t.t(end,1), unit_t);

time_tol = Time(0.01, 'ms');

if abs(end_time - topt.compute.tf) > time_tol
    disp(topt.compute)
    
    if runopt.draw_gui == true
        if mechguiopt.aborted == true
            dbclear if error
            error('mech_gui:abort', 'aborted by user')
        end
    end

    error_identifier = 'oc_mna:mainMNA:integrationError';
    error_str = 'Expected compute end time not reached. Probably due to an integration error.';
    
    fileID = fopen(fullfile(runopt.path.oc_mna, 'FAILURE'), 'w');
    fprintf(fileID, [error_str, '\n']);
    fprintf(fileID, 'Requested time: %s\n', end_time);
    fprintf(fileID, 'Computed time: %s\n', topt.compute.tf);
    fclose(fileID);
    
    error(error_identifier, error_str)
end

% save( fullfile( runopt.path.oc_mna, 'MNAparameters.mat'), ...
%     'mnaopt', '-v7.3');
%
% mnaParametersFile = matfile(fullfile( runopt.path.oc_mna, 'MNAparameters.mat'));

vohc_ss = mechopt.vohc_ss;
vihc_ss = mechopt.vihc_ss;

fn = sprintf('ss_voltage/%s.mat', runopt.hash.oc_mna);
dbprintf('Saving ss to %s\n', fn)
save(fn, 'vihc_ss', 'vohc_ss', '-v7.3')


%% Move data from scratch

files = {'Y.mat', 't.mat', 'cp_Y.mat', 'cp_t.mat'};

if use_scratch
    timer_copy = tic();
    if ~exist(runopt.path.oc_mna, 'dir')
        mkdir(runopt.path.oc_mna)
    end
    for j = 1:numel(files)
        if exist(fullfile(DIR, files{j}), 'file')
            % Currently, there are some problems with movefile in linu.
%             [status,message,messageId] = movefile( ...
%                 fullfile(DIR, files{j}), ...
%                 fullfile(runopt.path.oc_mna, files{j}));
%             if status == false
%                 error('Copy of file from scrath was not successfull.\nsource: %s\ndestination: %s\nmessage: %s\nmessage id: %s\n', ...
%                     fullfile(DIR, files{j}), ...
%                     fullfile(runopt.path.oc_mna, files{j}), ...
%                     message, messageId)
%             end
            src_f = fullfile(DIR, files{j});
            dst_f = fullfile(runopt.path.oc_mna, files{j});
            
            assert(~contains(src_f, '"'))
            assert(~contains(dst_f, '"'))
            
            cmd = sprintf('cp "%s" "%s"', src_f, dst_f);

            [status, message] = system(cmd);

            delete(fullfile(DIR, files{j}));
            if status ~= 0
                error('Copy of file from scrath was not successfull.\nsource: %s\ndestination: %s\nmessage: %s\n', ...
                    fullfile(DIR, files{j}), ...
                    fullfile(runopt.path.oc_mna, files{j}), ...
                    message)
            end
        else
            warning('[Copying files from scratch]: The file %s does not exist!', fullfile(DIR, files{j}))
        end
    end
    dbprintf('Copy time from scratch: %s', disp_toc(toc(timer_copy)));
end

%%

fclose('all');

function pos_ind = circuit_pos_ind(mnaopt)
    
    K = mnaopt.circuit.num_v + mnaopt.circuit.num_node;

    pos_ind = zeros(1,K*mnaopt.Numstacks);
    
    for i = 1:mnaopt.circuit.num_node
        pos_ind(i:K:end) = 1:mnaopt.Numstacks;
    end
    for i = (mnaopt.circuit.num_node+1):K
        pos_ind(i:K:end) = 1:mnaopt.Numstacks;
    end

end

function abs_tol = circuit_err_tol(mnaopt, stimulus, loc_scale)
    
    K = mnaopt.circuit.num_v + mnaopt.circuit.num_node;

    abs_tol = zeros(1,K*mnaopt.Numstacks);
    
    for i = 1:mnaopt.circuit.num_node
        abs_tol(i:K:end) = mnaopt.volt_tol.(mnaopt.simulation_units.voltage) .* loc_scale;
    end
    for i = (mnaopt.circuit.num_node+1):K
        % tolerance for single IHC or 3 OHCs
        abs_tol(i:K:end) = mnaopt.curr_tol.(mnaopt.simulation_units.current) * mnaopt.NumstacksScaling .* loc_scale;
    end
    
    % scale with db SPL    
    % db_scale = max([1, stimulus.amplitude]) / 10;
    %
    % if stimulus.amplitude is lower than ref_spl, then db_scale is < 1
    %                    higher                               > 1
    ref_spl = mnaopt.tol_ref_spl;
    db_scale = 10^((max([stimulus.amplitude]) - ref_spl)/20);
    
    abs_tol = abs_tol * db_scale;

    dbprintf('voltage abstol = %g V\n', mnaopt.volt_tol.(mnaopt.simulation_units.voltage) * db_scale)
    dbprintf('current abstol = %g A\n', mnaopt.curr_tol.(mnaopt.simulation_units.current) * db_scale)

end

function abs_tol = channel_err_tol(mnaopt, stimulus, loc_scale)
    
    % char_size = [
    %    1e-4, ... p
    %     6e-4, ... dp/dt
    % ];

    char_size = [
        0.1, ... p
        600, ... dp/dt
    ];

    tol0 = [
        1e-3, ... p
        1e-2, ... dp/dt
    ];
    
    % scale with db SPL    
    % db_scale = max([1, stimulus.amplitude]) / 10;
    %
    % if stimulus.amplitude is lower than ref_spl, then db_scale is < 1
    %                    higher                               > 1
    % ref_spl = 40;
    % ref_spl = 0;
    ref_spl = 60;
    db_scale = 10^((max([stimulus.amplitude]) - ref_spl)/20);
    
    char_size = tol0 .* char_size * db_scale;
    
    assert(numel(char_size) == 2)

    dbprintf('popen abstol = %g\n', char_size(1))
    dbprintf('dpopen/dt abstol = %g\n', char_size(2))

    % scale for velocity
    % [char_frequency, ~] = frequency_map();
    [~, char_frequency] = characteristic_frequency_model('Li');
    velocity_scale = char_frequency(mnaopt.xgrid);
    ref_freq = 1000; % Hz
    velocity_scale = max(1, velocity_scale/ref_freq);
    % velocity_scale = ones(1, mnaopt.Numstacks);
    
    % channels
    tol_channel_prob = char_size(1) * loc_scale;
    tol_channel_prob_der = char_size(2) * velocity_scale .* loc_scale;
    
    abs_tol = [];
    
    for i = 1:numel(mnaopt.channels)
        tol = [tol_channel_prob; tol_channel_prob_der];
        abs_tol = [abs_tol, tol(:)'];
    end
end

function pos_ind = channel_pos_ind(mnaopt)
    
    pos_ind = [];
    
    for i = 1:numel(mnaopt.channels)
        tmp = [1:mnaopt.Numstacks; 1:mnaopt.Numstacks];
        pos_ind = [pos_ind, tmp(:)'];
    end
end

function abs_tol = dvohc_err_tol(mnaopt, stimulus, loc_scale)
    ref_spl = 0;
    db_scale = 10^((max([stimulus.amplitude]) - ref_spl)/20);

    % [char_frequency, ~] = frequency_map();
    [~, char_frequency] = characteristic_frequency_model('Li');
    velocity_scale = char_frequency(mnaopt.xgrid);
    ref_freq = 1000; % Hz
    velocity_scale = max(1, velocity_scale/ref_freq);
    % velocity_scale = ones(1, mechopt.Numstacks);

    char_size = 1e-5;
    char_size = char_size * db_scale;

    dbprintf('dvohc abstol = %g\n', char_size(1))

    abs_tol = char_size * velocity_scale .* loc_scale;

end

function pos_ind = dvohc_pos_ind(mnaopt)

    pos_ind = 1:mnaopt.Numstacks;

end

function abs_tol = mech_err_tol(mechopt, stimulus, loc_scale)
    
    % in nm
    char_size = [
        1, ... bmx
        6, ... bmv
        1 / 25, ... tmx
        6 / 25, ... tmv
    ];    

    char_size = [
        0.1, ... bmx
        1, ... bmv
        4*0.01, ... tmx
        4*0.1, ... tmv
    ];

    tol_factor = 1;

    char_size = char_size * tol_factor;

    % unscaled velocity ~ 6x higher than displacement
    % TM ~ 25x lower than BM
    
    % scale with db SPL    
    % db_scale = max([1, stimulus.amplitude]) / 10;
    %
    % if stimulus.amplitude is lower than ref_spl, then db_scale is < 1
    %                    higher                               > 1
    % ref_spl = 100;
    ref_spl = 80;
    % ref_spl = 0;
    db_scale = 10^((max([stimulus.amplitude]) - ref_spl)/20);
    
    char_size = char_size * db_scale;
    
    dbprintf('bmx abstol = %g nm\n', char_size(1))
    dbprintf('tmx abstol = %g nm\n', char_size(3))
    
    % in model units
%     char_size = char_size ./ [ ...
%         mechopt.NMkonst_BM, ...
%         mechopt.NMkonst_OHC_cilia, ...
%         mechopt.NMkonst_BM, ...
%         mechopt.NMkonst_OHC_cilia, ...
%     ];
    
    % scale velocity with frequency
    % [char_frequency, ~] = frequency_map();
    [~, char_frequency] = characteristic_frequency_model('Li');
    velocity_scale = char_frequency(mechopt.xgrid);
    ref_freq = 1000; % Hz
    velocity_scale = max(1, velocity_scale/ref_freq);
    % velocity_scale = ones(1, mechopt.Numstacks);

    abs_tol = [ ...
        char_size(1) * loc_scale, ...
        char_size(2) * velocity_scale .* loc_scale, ...
        char_size(3) * loc_scale, ...
        char_size(4) * velocity_scale .* loc_scale];
end

function pos_ind = mech_pos_ind(mechopt)
    
    pos_ind = [ ...
        1:mechopt.Numstacks, ...
        1:mechopt.Numstacks, ...
        1:mechopt.Numstacks, ...
        1:mechopt.Numstacks];
end

end
