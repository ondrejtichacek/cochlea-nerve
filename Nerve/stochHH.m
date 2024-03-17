function [y, t, i, n, Na_max, calcN_Na, mh, initN_Na] = stochHH(modelname, t_input, I_input, replica, topt, hhopt, runopt, memopt)

% stochasticHH runs a simulation of a stochastic Hodgkin-Huxley model of:
%
%    Bruce, I. C. (2007). "Implementation issues in approximate methods for
%    stochastic Hodgkin???Huxley models," Annals of Biomedical Engineering
%    35(2):315-318.
%
% [V,N_Na,Istim,t] = stochasticHH(modelname,t_input,I_input)
%
% V is the relative transmembrane potential in units of millivolt.
% N_Na is the number of open sodium channels.
% Istim is the stimulus waveform in units of microampere.
% t is the time in units of millisecond.
%
% modelname is the model type, which may be one of the following:-
%    - 'hhfloor': Deterministic Hodgkin-Huxley model with N_Na rounded down.
%    - 'hhnint':  Deterministic Hodgkin-Huxley model with N_Na rounded to
%                 nearest integer.
%    - 'fxfloor': Fox algorithm with N_Na rounded down.
%    - 'fxnint':  Fox algorithm with N_Na rounded to nearest integer.
%    - 'cw':      Chow & White algorithm.
% I0 is the current amplitude in units of picoampere.  A single value
%    (e.g., I0 = 22) produces a single pulse of duration 100 microsecond
%    and amplitude I0 pA.
%    A pair of values (e.g., I0 = [22 35]) produces a conditioning prepulse
%    of duration 500 microsecond & amplitude I0(1) pA and a pulse of
%    duration 100 microsecond & amplitude I0(2) pA.
%
% Ian C. Bruce (ibruce@ieee.org), Faheem Dinath & Melissa T. Perri (c) 2007

%%

% Constant variables used in simulations. Corrected values are given in
% Table 1 of Bruce (2007)

gNa = 2.569e-8;     % maximum sodium conductance in mS [gamma_Na = 25.69 pS/channel]
Na_max = 1000;      % total number of sodium channels
Cm = 0.0714e-6;     % membrane capacitance in uF [Cm = 0.0714nF]
Rm = 1953.49e3;     % membrane resistance in kOhms [Rm = 1.95349e3 MOhms]
V_rest = -78;       % membrane resting potential in mV
E_Na = 66 - V_rest; % relative sodium Nernst potential in mV

%% Convrert units to doubles

unit = 'ms';

%% Configure solver and problem functions

dt = 1 / hhopt.samplingFrequency;
dt = dt.(unit);

%%

% Set up the appropriate functions for the chosen model
switch modelname
%     case 'hhfloor'
%         calcN_Na = @hhfloor;
%         initN_Na = @inithhfloor;
    case 'hhnint'
        mh = @hh_mh;
        calcN_Na = @hhnint;
        initN_Na = @hh_mh_init;
%     case 'fxfloor'
%         calcN_Na = @fxfloor;
%         initN_Na = @initfxfloor;
%         global m h;
%     case 'fxnint'
%         calcN_Na = @fxnint;
%         initN_Na = @initfxnint;
%         global m h;
    case 'cw'
        hs = hstruct('Nmh', []);
        hs = nmh_struct();
        mh = @hh_mh_mock;
        calcN_Na = @(m, h, markovrates, Na_max) cw(m, h, markovrates, Na_max, dt, hs);
        initN_Na = @(markovrates) initcw(markovrates, Na_max, dt, hs);
        % global Nmh;
    otherwise
        error('Model type "%s" unknown', modelname)
end

%% Init result files

resPath = replica.dir;

resFiles = {'i','y','t'};

switch hhopt.save_method
    case 'matlab_matfile'
        ext = 'mat';
        
    case {'matlab_fwrite', 'c_posix'}                
        ext = 'dat';
end

switch hhopt.save_method
    case 'struct'
        tmp = {};
        for i = 1:numel(resFiles)
            tmp{2*i-1} = resFiles{i};
            tmp{2*i} = struct;
        end

        res = hstruct(tmp{:});

        saving_to_file = false;
    otherwise

        saving_to_file = true;

        for j = 1:length(resFiles)
            r = resFiles{j};
            res.(r) = WrapFile( ...
                'name', r, ...
                'extension', ext, ...
                'dir', resPath, ...
                'type', 'double', ...
                'save_method', hhopt.save_method ...
                );
        
            checkpoints.(r) = WrapFile( ...
                'name', ['cp_', r], ...
                'extension', ext, ...
                'dir', resPath, ...
                'type', 'double', ...
                'save_method', hhopt.save_method ...
                );
            
            if strcmp(replica.action, 'create')
                if exist(res.(r).path, 'file')
                    warning('Deleting old results file: %s', res.(r).path);
                    sdelete(res.(r).path);
                end
                if exist(checkpoints.(r).path, 'file')
                    warning('Deleting old checkpoints file: %s', checkpoints.(r).path);
                    sdelete(checkpoints.(r).path);
                end
            end
        end
end

%% Configure solver and problem functions

% tspan = [t_input(1), t_input(end)];
tspan = t_input;

SOLVER = @odeEuler;
ODEFUN = @derHH;

solveropt = solverOpt( ...
    'Precompute', @(tt) interp1( t_input.(unit), I_input, tt ), ...
    'TimeStep', dt ...
    );

if runopt.verbose >= 3
    ws = WaitbarStorage();
    ws.num_div = 1;
    
    if runopt.waitbarFunctionAvailable
        tf = tspan(end).(unit); % faster to take out of anon. fcn
        OutputFcn = @(t, y, flag) odewaitbar(t, y, flag, tf, ws);
    else
        tf = tspan(end).(unit); % faster to take out of anon. fcn
        OutputFcn = @(t, y, flag) odewaittext(t, y, flag, tf, ws);
    end
    solveropt.OutputFcn = OutputFcn;
    solveropt.OutputFcnEvalInterval = 10*dt;
end

%% Open result files for writing

switch hhopt.save_method
    case 'matlab_matfile'    
        res.t = matfile(res.t.path, 'Writable', true);
        res.y = matfile(res.y.path, 'Writable', true);
        checkpoints.t = matfile(checkpoints.t.path, 'Writable', true);
        checkpoints.y = matfile(checkpoints.y.path, 'Writable', true);
        
    case 'matlab_fwrite'
        res.t.fopen('w');
        res.y.fopen('w');
        checkpoints.t.fopen('w');
        checkpoints.y.fopen('w');

    case 'struct'
        checkpoints.t = struct();
        checkpoints.y = struct();
end

%% Initial conditions

% n ... number of slices
n = size(I_input,2);

switch replica.action
    case 'create'
        V0 = zeros(n,1);

        % Calculate sodium gating particle transition rates for initial (resting)
        % membrane potential
        markovrates0 = calcrates(V0);

        % Calculate initial ???
        [m0, h0] = initN_Na(markovrates0);

        y0 = [ ...
            V0; ... membrane potential
            m0; ... fraction of sodium activation particles
            h0 ... fraction of sodium inactivation particles
            ];

    case 'extend'
        y0 = res.y.y(end,:);
        y0 = y0';
end

tspan = tspan.(unit);

%% Solve problem

disp_title(sprintf('Solving HH'), '-');

wrapopt.append = 'true';
wrapopt.res = res;
wrapopt.verbose = runopt.verbose;
wrapopt.NumDiv = hhopt.NumDiv;
wrapopt.memopt = memopt;
wrapopt.save_method = hhopt.save_method;
wrapopt.save_checkpoints = true;
wrapopt.checkpoint_period = topt.compute.checkpoint_period.(unit);
wrapopt.checkpoints = checkpoints;

% Get the length of saved checkpoints
if saving_to_file && exist(checkpoints.t.Properties.Source, 'file')
    wrapopt.checkpoints_last_t_index = max(matfileWhosByName(checkpoints.t, 't', 'size')); % = memmory friendly length(checkpoints.t.t)
else
    wrapopt.checkpoints_last_t_index = 0;
end

[ statistics ] = odeWrapper( wrapopt, SOLVER, ODEFUN, tspan, y0, solveropt, ...
    I_input, t_input.(unit), Na_max, Rm, gNa, E_Na, Cm, calcN_Na, mh, dt );


res.t.unit = unit;

%% Close result files

if strcmp(hhopt.save_method, 'matlab_fwrite')
    res.t.fclose();
    res.y.fclose();
end

%% Write signal and finalize

T = Time(res.t.t, res.t.unit);
T = T(T >= t_input(1) & T <= t_input(end) + Time(1000*eps, 'ms'));

nX = size(I_input, 2);

switch hhopt.save_method
    case 'matlab_matfile'
        res.i = matfile(res.i.path, 'Writable', true);
        
        if exist(res.i.Properties.Source, 'file')
            n0 = matfileWhosByName(res.i, 'Istim', 'size');
            
            n0t = n0(1);
            n0t = n0t - 1;
            
            n0x = n0(2);
            assert(nX == n0x, 'Spatial dimension missmatch!')
        else
            n0t = 0;
        end
        
        tind = n0t + (1:length(T));
        
        if nX == 1
            % res.i.Istim(tind, 1:nX) = interp1(t_input.us, I_input, T.us); % Stimulus vector
            F = griddedInterpolant(t_input.us, I_input);
            res.i.Istim(tind, 1) = F(T.us);
        else
            F = griddedInterpolant({t_input.us, 1:nX}, I_input);
            res.i.Istim(tind, 1:nX) = F({T.us, 1:nX});
        end
        
        clear t_input I_input
    case 'struct'
        
        n0t = 0;
        
        tind = n0t + (1:length(T));
        
        if nX == 1
            % res.i.Istim(tind, 1:nX) = interp1(t_input.us, I_input, T.us); % Stimulus vector
            F = griddedInterpolant(t_input.us, I_input);
            res.i.Istim(tind, 1) = F(T.us);
        else
            F = griddedInterpolant({t_input.us, 1:nX}, I_input);
            res.i.Istim(tind, 1:nX) = F({T.us, 1:nX});
        end
        
        clear t_input I_input

    case {'matlab_fwrite', 'c_posix'}
        
        error('Not implemented');
end

y = res.y;
t = res.t;
i = res.i;

%% Gather results

% V = y(:,1:n); % Relative transmembrane potential vector

% m = y(:,n+(1:n)); % Fraction of sodium activation particles
% h = y(:,2*n+(1:n)); % Fraction of sodium inactivation particles

% N_Na = calcN_Na(m, h, Na_max); % Open sodium channel count vector

% Istim = interp1(t_input,I_input,t); % Stimulus vector
