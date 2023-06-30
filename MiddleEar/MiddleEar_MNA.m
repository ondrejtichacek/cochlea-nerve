function [t, y, y_desc] = MiddleEar_MNA(topt, stimulus, midopt, runopt, opt)
arguments
    topt     (1,1) timeOpt
    stimulus (1,1) {mustBeA(stimulus, ["Signal", "CompoundSignal"])}
    midopt   (1,1) midOpt
    runopt   (1,1) runOpt
    opt      (1,1) globalOpt
end


circuit = midopt.circuit;
simulation_units = midopt.simulation_units;

%%

dB = stimulus.amplitude;

p0 = 20e-6; % 20 uPa
p = p0 * 10.^(dB/20);

AM = p * midopt.pv_TM; % scaling factor at tympanic membrane

%% Input signal

if numel(stimulus) == 1
    signal_fun = stimulus.eval_fcn(simulation_units.time);
else
    error('Not fully tested')
    for i = 1:numel(stimulus)
        signal_fun_partial{i} = stimulus(i).eval_fcn(simulation_units.time);
    end 
    signal_fun = @(t) sum(cellfun(@(f) f(t), signal_fun_partial, ...
        'UniformOutput', true));
end

assert(all(size(AM) == size(signal_fun(0))))

% Not much faster
%
% frequency_unit = Time.inverse_unit(simulation_units.time);
% omega = stimulus.om.(frequency_unit);
% envelope = stimulus.get_envelope_fun(simulation_units.time);
% 
% signal_fun_fast = @(t) envelope(t) .* sin(omega .* t);
% 
% t = tic;
% for tt = linspace(tspan(1), tspan(2), 1000000)
%     signal_fun(tt);
% end
% toc(t)
% t = tic;
% for tt = linspace(tspan(1), tspan(2), 1000000)
%     signal_fun_fast(tt);
% end
% toc(t)

%% Assemble and solve MNA equations

[ A, C, Z, y_desc ] = RLC(...
    circuit.resistors, ...
    circuit.inductors, ...
    circuit.capacitors, ...
    circuit.vsources, ...
    circuit.isources, ...
    circuit.num_node);

% Y0 = zeros(size(Z));
% Y0 = linsolve(A, Z, struct('RECT', true));
Y0 = lsqminnorm(A, Z);
% A*Y0 - Z

Z0 = Z;
Z0(end) = 0;

Z_signal = zeros(size(Z));
Z_signal(end) = Z(end);

JACOBFUN = @(t,y) -A;
% ODEFUN_INIT = @(t,y) -A * y + Z*0;
ODEFUN_INIT = @(t,y) -A * y + Z0;
ODEFUN = @(t,y) -A * y + Z0 + Z_signal * sum(AM .* signal_fun(t));

t0 = topt.total.t0;
tspan = [t0.(simulation_units.time),topt.total.tf.(simulation_units.time)];

reltol = midopt.ode_reltol;

%% Configure solver

if isa(stimulus, 'PureTone')
    maxstep = 1/max([stimulus.frequency]) / 32;
elseif isa(stimulus, 'Click') || isa(stimulus, 'Click2')
    maxstep = midopt.ode_maxstep;
else
    maxstep = midopt.ode_maxstep;
end

initialstep = midopt.ode_initialstep;
initialstep = min([maxstep, initialstep]);

% Absolute tolerance
abs_tol = circuit_err_tol(midopt, stimulus);

solveropt = odeset( ...
    'Mass', C, ...
    'MassSingular', 'yes', ...
    'MStateDependence', 'none', ...
    'Jacobian', JACOBFUN, ...
    'MaxStep', maxstep.(simulation_units.time), ...
    'InitialStep', initialstep.(simulation_units.time), ...
    'AbsTol', abs_tol, ...
    'RelTol', reltol, ...
    'BDF', 'on', ...
    'MaxOrder', 2, ...
    ...'Vectorized', 'on', ...
    'Stats','on' ...
    );

%% Compute consistent initial conditions

[Y0_consistent, Yp0_consistent] = me_mna.InitialConditions( ...
    ODEFUN_INIT, solveropt, midopt, runopt, opt, [], Y0);

y0 = Y0_consistent;

%% Set tolerance according to envelope

tolerance_rel_envelope = false; % currently not working properly

if tolerance_rel_envelope == true

    F = envelope_tolerance(stimulus, tspan, ...
        simulation_units.time, ...
        stimulus.samplingFrequency, ...
        'amplifier', false, ...
        'plotflag', false);

    abs_tol = abs_tol / 10;

    abs_tol = abs_tol(:); % IMPORTANT

    abs_tol = @(t) F(t) * abs_tol;
end

solveropt = odeset(solveropt, 'AbsTol', abs_tol);

%%

tt = tic;
[t, y] = mtode.ode15s(ODEFUN, tspan, y0, ...
        odeset(solveropt, 'OutputFcn', [], 'Stats', 'on'));
disp_toc(toc(tt));

zero_ind = find(t >= 0, 1);
t = t(zero_ind:end, :);
y = y(zero_ind:end, :);


%%

function abs_tol = circuit_err_tol(midopt, stimulus)
    
    num_v = midopt.circuit.num_v;
    
    num_node = midopt.circuit.num_node;

    num_ind = midopt.circuit.num_ind;

    abs_tol = ones(num_node + num_ind + num_v, 1);
    
    ref_scale = [ ...
        1e-3; % VN_1
        1e-3; % VN_2
        1e-3; % VN_3
        1e-3; % VN_4
        1e-3; % VN_5
        1e-3; % VN_6
        1e-3; % VN_7
        1e-4; % VN_8
        1e-3; % VN_9
        1e-4; % VN_10
        1e-4; % VN_11
        1e-5; % VN_12
        ...
        1e-5; % VN_13
        1e-5; % VN_14
        1e-5; % VN_15
        1e-5; % VN_16
        1e-6; % VN_17
        1e-6; % VN_18
        ...
        1e-7; % IL_1
        1e-7; % IL_2
        1e-7; % IL_3
        1e-7; % IL_4
        1e-7; % IL_5
        1e-7; % IL_6
        ...
        1e-5; % IV_1
    ];

    % ref_scale = ones(size(ref_scale));

    abs_tol = 1e-6 .* abs_tol .* ref_scale;
    
    % scale with db SPL    
    % db_scale = max([1, stimulus.amplitude]) / 10;
    %
    % if stimulus.amplitude is lower than ref_spl, then db_scale is < 1
    %                    higher                               > 1

    ref_spl = 0;
    db_scale = 10^((max([stimulus.amplitude]) - ref_spl)/20);
    
    abs_tol = abs_tol * db_scale;

end

end