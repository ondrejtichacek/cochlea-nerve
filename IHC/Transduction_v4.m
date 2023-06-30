function [ y_out, t_out, V_out, y_dyn_out, n, tropt, f_ChannelsOpenSS, f_CalciumCurrent, f_TransmitterRelease ] = Transduction_v4( ...
        t, IHCVoltage, SR, replica, topt, antopt, runopt, memopt )
%TRANSDUCTION
%
%   t          ... vector of time points
%   IHCVoltage ... vector of voltage at t
%   SR         ... spontaneous rate => fibre type (string)
%   antopt     ... structure with options
%   runopt     ...
%   memompt    ...
%
%   S       ... structure with fields:
%     .T    ... vector of time points
%     .V    ... vector of voltage at T
%     .m    ... fraction of open Ca channels at T
%     .C    ... Ca concentration at T
%     .q    ... immediate store at T
%     .c    ... cleft at T
%     .w    ... reprocessing store at T
%     .m_inf
%     .I
%     .k
%     .parameters
%
% NOTE! t and T has generally different points but same range
%

%% INIT I/O FILES

resPath = replica.dir;
resFiles = {'V', 'y', 't', 'y_dyn'};
resFiles_IC = {'V_IC', 'y_IC', 't_IC'};

switch antopt.save_method
    case 'matlab_matfile'
        ext = 'mat';
    case {'matlab_fwrite', 'c_posix'}                
        ext = 'dat';        
end

[res, checkpoints] = prepare_res_files(replica, resFiles, ext, resPath, antopt);
[res_IC, checkpoints_IC] = prepare_res_files(replica, resFiles, ext, resPath, antopt);

%% SET UP OVERALL PARAMETERS

x_pos = antopt.positionsToExcite;

if numel(x_pos) ~= 1
    error("We currently support only one section at a time");
end

tropt = copy(antopt.transductionopt);

simulation_units = struct( ...
    'Voltage', 'V', ...
    'Frequency', 'Hz', ...
    'Concentration', 'M', ...
    'Time', 's');

dt = 1 / antopt.samplingFrequency;
dt = dt.(simulation_units.Time);
tspan = [t(1), t(end)];

channels_open_ss_parameters_normal = tropt.generate_channels_open_ss_parameters('normal', x_pos);
channels_open_ss_parameters_burst = tropt.generate_channels_open_ss_parameters('burst', x_pos);

transmitter_release_parameters = tropt.generate_transmitter_release_parameters(x_pos);

%% Background Ca_2+ concentration

% DOI: 10.1016/0196-0709(89)90030-6
% Ca concentration in endolymph
% C0 = 28.5e-6; % mean

C_Ca_background_base = 17.1e-6; % M
C_Ca_background_apex = 40.6e-6; % M

% interpolated
C_Ca_background = C_Ca_background_base + (C_Ca_background_apex - C_Ca_background_base) * x_pos;

C_Ca_background = C_Ca_background * tropt.C_Ca_background_factor;

%% CONFIGURE SOLVER AN PROBLEM FUNCTIONS

SOLVER = @odeEuler;
solveropt = solverOpt( ...
    'Precompute', @(tt) interp1( t.(simulation_units.Time), IHCVoltage, tt ), ...
    'TimeStep', dt ...
    );

if runopt.verbose >= 3
    ws = WaitbarStorage();
    ws.num_div = 1;

    if runopt.waitbarFunctionAvailable
        tf = tspan(end).(simulation_units.Time); % faster to take out of anon. fcn
        OutputFcn = @(t, y, flag) odewaitbar(t, y, flag, tf, ws);
    else
        tf = tspan(end).(simulation_units.Time); % faster to take out of anon. fcn
        OutputFcn = @(t, y, flag) odewaittext(t, y, flag, tf, ws);
    end
    solveropt.OutputFcn = OutputFcn;
    solveropt.OutputFcnEvalInterval = 10*dt;
end

tspan = tspan.(simulation_units.Time); % dropping the units

[numSteps, numSamples] = odeEuler_tspan(tspan, dt);

switch antopt.ant
    case 'v4'
        
        do_equlibration = false;
        
        % number of Ca_V1.3 channels in the vicinity of the synapse
        num_CaV13 = tropt.num_CaV13;
        
        tau_CaV13 = tropt.tau_CaV13.s;
        tau_CaV13_blocked =  tropt.tau_CaV13_blocked.s;
        
        num_vesicles = tropt.num_release_sites - tropt.num_inactive_release_sites;
        
        tmp = namedargs2cell(tropt.ribbon_synapse_properties);
        
        rs = RibbonSynapse_v4( tmp{:}, ...
            'plotflag', true, ...
            'num_release_sites', tropt.num_release_sites, ...
            'num_channels', tropt.num_CaV13 ...
                );
        
        % exportgraphics(gcf,'H19a.pdf','ContentType','vector')
        
        % inactivate some release sites
        if tropt.num_inactive_release_sites > 0

            active = sort(randi_unique_vec(tropt.num_release_sites, tropt.num_release_sites - tropt.num_inactive_release_sites));

            assert(numel(active) == numel(unique(active)))

            rs.xv = rs.xv(active,:);
            % rs.xc = rs.xc; % no change
            % rs.psi = rs.psi; % no change
            rs.rho = rs.rho(active,:);

        end

        channels = Channels(num_CaV13, tau_CaV13, tau_CaV13_blocked);

        % inital_p_open = 0.2;
        inital_p_open = 1;

        for iii = 1:channels.num
            if rand(1) < inital_p_open
                channels.mode(iii) = 1;
            else
                channels.mode(iii) = 0;
            end
        end

        % p_inactivated = 0.7;
        p_inactivated = 0.0;

        for iii = 1:channels.num
            if rand(1) < p_inactivated
                channels.state(iii) = 'i';
            end
        end

        V0t = channels_open_ss_parameters_normal{1};
        S0t = channels_open_ss_parameters_normal{2};

        channels.alpha = tropt.channels_markov_properties.alpha;
        channels.kp0 = tropt.channels_markov_properties.kp0;
        channels.V0t = V0t.V;
        channels.S0t = S0t.V;

        vesicles = Vesicles(num_vesicles);
        
        d = 3; %           ...  geometry factor
        D = tropt.Ca_diffusion_coefficient;
        
        % r = [tropt.d_channel_vesicle;  % distance vesicle--channel
        %      tropt.d_channel];        % distance of channel char. point from membrane
        
        % Relative tolerance for Ca_conc calculation. The lower the
        % rel_tol is, the longer history we need to calculate. Usually used
        % value is 0.05.
        rel_tol = tropt.Ca_conc_rel_tol;


        % Alternatively to rel_tol, we could set max_history
        % max_history = 10000;

        for i = 1:num_CaV13
            
            rho = rs.rho(:,i) * 1e-9; % channel mouth--vesicle membrane
            psi = rs.psi(:,i) * 1e-9; % channel mouth--channel mouth

            psi_nernst = sqrt(psi.^2 + tropt.d_nernst.^2); % channel mouth--channel nernst point (above channel mouth)
            
            % r = [tropt.d_nernst; rho];
            r = [psi_nernst; rho];

            ps(i) = Diffusion.PointSource(d, D, r, dt, 0:numSamples, ...
                'rel_tol', rel_tol);
                % 'max_history', max_history);
        end

        N_history = [ps.N];
        fprintf('History of avg. %d (min %d, max %d) steps required for %g rel error.\n', ...
            round(mean(N_history)), min(N_history), max(N_history), rel_tol)
        
        % assign channels to vesicles
        % by default assign all channels to all vesicles (correct, but comp. intesive)
        for i = 1:num_vesicles
            vesicles.close_channels{i} = 1:num_CaV13;
        end

        % alternatively we could assign only the close channels (not currently used1)
        calculate_only_nearest_channels = false;
        if calculate_only_nearest_channels
            num_channels_per_vesicle = num_CaV13 / num_vesicles;
            nc = 0;
            for i = 1:num_vesicles
                vesicles.close_channels{i} = nc + (1:num_channels_per_vesicle);
                nc = nc + num_channels_per_vesicle;
            end
        end

        % ODEFUN = @TransductionRHS_v4;
        ODEFUN = @TransductionRHS_v5;
        
    otherwise
        error('Unrecognized value %s', antopt.ant)
end

%% INITIAL CONDITIONS

% assumptions at t = 0:
%   - I0 = 0
%   - C0 = C_background


% n ... number of slices
n = size(IHCVoltage,2);

assert(n == 1); % for now

V_steady_state = IHCVoltage(1,:)';

% creating a new replica or extending an old one
switch replica.action
    case 'create'   % initial conditions for the new replica
        % V0 : initial voltage [V]
        %
        % m0 : Resting P_open of Cachannels - Boltzmann function
        %
        % C0 : initial Ca concentration
        %
        % q0 : free neurotransmitter (NT) pool - initially all NT is here
        % c0 : cleft - initially empty
        % w0 : reprocessing storage - initially empty
        

        % V0 = IHCVoltage(1,:)';
        % m0 = ChannelsOpenSS( Voltage(V0, 'V'), channels_open_ss_parameters_normal{:} );

        if isempty(antopt.initialConditions)

            % m0 = sum([channels.mode] == 1) / channels.num;
            m0 = sum([channels.state] == 'o') / channels.num;

            C0 = C_Ca_background;
            C0 = repmat(C0, 1, num_vesicles); % row
            
            q0 = ones(n, num_vesicles);  % all NT in free pool
            
            w0 = zeros(n,1);           % empty
        else

            [m0_mean, C0, q0, w0] = prepare_IC(antopt.initialConditions);
            [m0_mean, C0, q0, w0] = prepare_IC_new(antopt.initialConditions_ss);

            for iii = 1:channels.num
                if rand(1) < m0_mean
                    channels.state(iii) = 'o';
                else
                    channels.state(iii) = 'c';
                end
            end
        
            m0 = sum([channels.state] == 'o') / channels.num;

        end
        
        m_block0 = sum([channels.mode] == 0) / channels.num;

        I0 = 0;
        
        c_glut0 = zeros(n,1);           % empty        
        c_prot0 = zeros(n,1);           % empty

        y0 = [m0, m_block0, I0, C0, q0, c_glut0, w0, c_prot0]';
        
        state_inactivated = zeros(n,1);
        state_normal = num_CaV13 * ones(n,1);
        state_burst = zeros(n,1);
        
        y0 = [y0; state_inactivated; state_normal; state_burst];
        
        size_info = create_size_info( ...
            "CaV13_channels_fraction", 1, ...
            "CaV13_channels_blocked", 1, ...
            "Ca_current", 1, ...
            "Ca_concentration", num_vesicles, ...
            "NT_free", num_vesicles, ...
            "NT_cleft", 1, ...
            "NT_reprocessing", 1, ...
            "proton_cleft", 1, ...
            "CaV13_num_inactivated", 1, ...
            "CaV13_num_normal", 1, ...
            "CaV13_num_burst", 1);

        size_info_2 = create_size_info( ...
            "NT_free", num_vesicles, ...
            "NT_cleft", 1, ...
            "NT_reprocessing", 1, ...
            "proton_cleft", 1);
        
    case 'extend'
        y0 = res.y.y(end,:);
        y0 = y0';
end

% parameters from Sumner 2002
if isscalar(tropt.x), tropt.x = repmat(tropt.x,[n,1]); end
if isscalar(tropt.y), tropt.y = repmat(tropt.y,[n,1]); end
if isscalar(tropt.r), tropt.r = repmat(tropt.r,[n,1]); end
if isscalar(tropt.l), tropt.l = repmat(tropt.l,[n,1]); end

uconv = @(var) Unit.batch_convert(simulation_units, var);

%% SOLVE PROBLEM

if runopt.verbose >= 3
    disp_title(sprintf('Solving ODE for %s', SR), '-');
end

wrapopt.append = 'true';
wrapopt.res = res;
wrapopt.verbose = runopt.verbose;
wrapopt.NumDiv = antopt.NumDiv;
wrapopt.memopt = memopt;
wrapopt.save_method = antopt.save_method;
wrapopt.save_checkpoints = true;
wrapopt.checkpoint_period = topt.compute.checkpoint_period.(simulation_units.Time);
wrapopt.checkpoints = checkpoints;

% Get the length of saved checkpoints
if exist(checkpoints.t.Properties.Source, 'file')
    wrapopt.checkpoints_last_t_index = max(matfileWhosByName(checkpoints.t, 't', 'size')); % = memmory friendly length(checkpoints.t.t)
else
    wrapopt.checkpoints_last_t_index = 0;
end

if strcmp(replica.action, 'create') && do_equlibration
    
    % Do a short equilibration
    
    tspan_IC = Time([0, 1], 'ms');
    tspan_IC = tspan_IC.(simulation_units.Time);
    
    solveropt_IC = solverOpt( ...
        'Precompute', @(tt) repmat(IHCVoltage(1,:), [numel(tt), 1]), ...
        'TimeStep', dt ...
        );
    
    wrapopt_IC = wrapopt;
    wrapopt_IC.res = res_IC;
    wrapopt_IC.checkpoints = checkpoints_IC;
    
    [ statistics_IC ] = odeWrapper( wrapopt_IC, SOLVER, ODEFUN, tspan_IC, y0, solveropt_IC, ...
        size_info, channels, vesicles, ps, tropt.y.Hz, tropt.l.Hz, tropt.x.Hz, tropt.r.Hz, ...
        tropt.G_Ca, C_Ca_background, ...
        uconv(channels_open_ss_parameters_normal), ...
        uconv(channels_open_ss_parameters_burst), ...
        uconv(transmitter_release_parameters), V_steady_state, dt );
   
    y0 = res_IC.y.y(end,:);
    y0 = y0';
    
end

% vesicles.reinit(num_vesicles, numSamples);

[ statistics ] = odeWrapper( wrapopt, SOLVER, ODEFUN, tspan, y0, solveropt, ...
    size_info, channels, vesicles, ps, tropt.y.Hz, tropt.l.Hz, tropt.x.Hz, tropt.r.Hz, ...
    tropt.G_Ca, C_Ca_background, ...
    uconv(channels_open_ss_parameters_normal), ...
    uconv(channels_open_ss_parameters_burst), ...
    uconv(transmitter_release_parameters), V_steady_state, dt );

res.t.unit = simulation_units.Time;

res.y.size_info = size_info;

%% GATHER RESULTS

T = Time(res.t.t, res.t.unit);
T = T(T >= t(1) & T <= t(end));

nX = size(IHCVoltage, 2);

%%

% debug_plot = true;
debug_plot = false;

if debug_plot
    plot_this()
end


%%

switch antopt.save_method
    case 'matlab_matfile'
        
        if exist(res.V.Properties.Source, 'file')
            n0 = matfileWhosByName(res.V, 'V', 'size');
            
            n0t = n0(1);
            n0t = n0t - 1;
            
            n0x = n0(2);
            assert(nX == n0x, 'Spatial dimension missmatch!')
        else
            n0t = 0;
        end
        
        tind = n0t + (1:length(T));
        
        if nX == 1
            % res.V.V(tind, 1:nX) = interp1(t.us, Voltage, T.us); % Stimulus vector
            F = griddedInterpolant(t.us, IHCVoltage);
            res.V.V(tind, 1) = F(T.us);
        else
            F = griddedInterpolant({t.us, 1:nX}, IHCVoltage);
            res.V.V(tind, 1:nX) = F({T.us, 1:nX});
        end
        
        clear Time Voltage
                
    case {'matlab_fwrite', 'c_posix'}
        
        error('Not implemented');
end

res.V.unit = simulation_units.Voltage;

res.y.vesicle_release_events = vesicles.release_events(~isnan(vesicles.release_events));

%% RETURN VALUES
y_out = res.y;
t_out = res.t;
V_out = res.V;

f_ChannelsOpenSS = @(V) CaV13.ChannelsOpenSS(V, channels_open_ss_parameters_normal{:});

% f_ChannelsOpenSS = @(V) CaV13.ChannelsOpenSS(V, channels_open_ss_parameters_normal{:});    

f_CalciumCurrent = @(V, m, E) CalciumCurrent(V.V, m, tropt.G_Ca, E.V);

tr_par = uconv(transmitter_release_parameters);
f_TransmitterRelease = @(C) RibbonSynapse_v4.TransmitterRelease(C, tr_par{:});

%%

C_vesicles = decompose_z_all_t(y_out.y, 'Ca_concentration', size_info);
tt = Time(t_out.t, t_out.unit);

F = griddedInterpolant(tt.(simulation_units.Time), C_vesicles);

solveropt = solverOpt( ...
    ...'Precompute', @(tt) interp1( t.(simulation_units.Time), C_vesicles, tt ), ...
    'Precompute', @(tt) F(tt), ...
    'TimeStep', dt ...
    );

ODEFUN = @NTdynamicsRHS_v5;

n_dyn_rep = 100;

y_dyn = cell(1, n_dyn_rep);
y_dyn_release = cell(1, n_dyn_rep);

for i = 1:n_dyn_rep

    if ~isempty(antopt.initialConditions) % random IC every time
        [~, ~, q0, w0] = prepare_IC(antopt.initialConditions);
        [~, ~, q0, w0] = prepare_IC_new(antopt.initialConditions_ss);
    end
    
    y0 = [q0, c_glut0, w0, c_prot0]';

    vesicles = Vesicles(num_vesicles);

    [t_dyn, y_dyn{i}] = odeEuler(ODEFUN, tspan, y0, solveropt, size_info_2, channels, vesicles, ps, tropt.y.Hz, tropt.l.Hz, tropt.x.Hz, tropt.r.Hz, ...
        tropt.G_Ca, C_Ca_background, ...
        uconv(channels_open_ss_parameters_normal), ...
        uconv(channels_open_ss_parameters_burst), ...
        uconv(transmitter_release_parameters), V_steady_state, dt );
    
    y_dyn_release{i} = vesicles.release_events(~isnan(vesicles.release_events));

end

res.y_dyn.y = y_dyn;
res.y_dyn.size_info = size_info_2;
res.y_dyn.vesicle_release_events = y_dyn_release;

y_dyn_out = res.y_dyn;

% [ statistics ] = odeWrapper( wrapopt, SOLVER, ODEFUN, tspan, y0, solveropt, ...
%     size_info, channels, vesicles, ps, tropt.y.Hz, tropt.l.Hz, tropt.x.Hz, tropt.r.Hz, ...
%     tropt.G_Ca, C_Ca_background, ...
%     uconv(channels_open_ss_parameters_normal), ...
%     uconv(channels_open_ss_parameters_burst), ...
%     uconv(transmitter_release_parameters), V_steady_state, dt );


%%

function [] = plot_this()
%TRANSDUCTION_V4_PLOT 

figure
c = zeros(numSamples+1, vesicles.num);

for ii = 1:vesicles.num
    if ii == 1
        subplot(3,1,1)
    elseif ii == 6
        subplot(3,1,2)
    elseif ii == 11
        subplot(3,1,3)
    end
    hold on
        
    for j = 1:channels.num
        c(:,ii) = c(:,ii) + ps(j).concentration(ii+1,:)';
    end
    plot(T.ms, c(1:end-1,ii)*1e6, 'DisplayName', sprintf('ves %d', ii));
    ylabel('Ca2+ concentration @ vesicle (uM)')
    xlabel('Time (ms)')
    legend('Location', 'NorthWest')
end

threshold = 250e-6;

figure
hold on
histogram(c*1e6, 'Normalization', 'probability')
histogram(c(c>threshold)*1e6, 'Normalization', 'probability')
xlabel('Ca2+ concentration (uM)')


tmp = {ps.concentration};
cc = cat(3,tmp{:});
cc = cc(2:end,:,:);
w = cc;
% w = w ./ sum(w(:));
[histw, ~, bin_centers] = histwv(cc(:), w(:), 'nbins', 256);

figure
hold on
plot(bin_centers*1e6, histw, ...
    'DisplayName', 'all events')
xlabel('Ca2+ concentration (uM)')

ww = zeros(size(cc));
for ii = 1:vesicles.num
    for j = 1:channels.num
        ww(ii,:,j) = c(:,ii) > threshold;
    end
end

w = cc.*ww;
% w = w ./ sum(w(:));
[histw, ~, bin_centers] = histwv(cc(:), w(:), 'nbins', 256);

plot(bin_centers*1e6, histw, ...
    'DisplayName', 'over threshold')
xlabel('Ca2+ concentration (uM)')

www = zeros(size(cc));
for ii = 1:vesicles.num
    for j = 1:channels.num
        www(ii,:,j) = c(:,ii) > threshold & all(cc(ii,:,:) < threshold, 3)';
    end
end

w = cc.*www;
% w = w ./ sum(w(:));
[histw, ~, bin_centers] = histwv(cc(:), w(:), 'nbins', 256);

plot(bin_centers*1e6, histw, ...
    'DisplayName', 'over threshold small contributions')
xlabel('Ca2+ concentration (uM)')

legend()

title('Ca2+ concentration @ vesicle')
subtitle('weighted by concentration value')

end

end

function [res, checkpoints] = prepare_res_files(replica, resFiles, ext, resPath, antopt)
for j = 1:length(resFiles)
    r = resFiles{j};
    res.(r) = WrapFile( ...
        'name', r, ...
        'extension', ext, ...
        'dir', resPath, ...
        'type', 'double', ...
        'save_method', antopt.save_method ...
        );

    checkpoints.(r) = WrapFile( ...
        'name', ['cp_', r], ...
        'extension', ext, ...
        'dir', resPath, ...
        'type', 'double', ...
        'save_method', antopt.save_method ...
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

    switch antopt.save_method
        case 'matlab_matfile'
            res.(r) = matfile(res.(r).path, 'Writable', true);
            checkpoints.(r) = matfile(checkpoints.(r).path, 'Writable', true);
        case 'matlab_fwrite'
            res.(r).fopen('w');
            checkpoints.(r).fopen('w');
    end
end

end

function v = decompose_z_all_t(z, variable, size_info)

si = size_info.(variable);
v = z(:, si.start : si.end);

end

function [m0, C0, q0, w0] = prepare_IC(IC_mean)
    %% Process IC (random based on mean)

    q0 = zeros(size(IC_mean.q0_mean));
    for i = 1:numel(q0)
        q0(i) = floor(IC_mean.q0_mean(i));
        q0(i) = q0(i) + (rand(1) < (IC_mean.q0_mean(i) - q0(i)));
    end

    w0 = zeros(size(IC_mean.w0_mean));
    for i = 1:numel(w0)
        w0(i) = floor(IC_mean.w0_mean(i));
        w0(i) = w0(i) + (rand(1) < (IC_mean.w0_mean(i) - w0(i)));
    end

    C0_mean = IC_mean.C0_mean;
    m0_mean = IC_mean.m0_mean;

    m0 = m0_mean;
    C0 = C0_mean;
end

function [m0, C0, q0, w0] = prepare_IC_new(IC_mean)
    %% Process IC (random based on random time)

    ind_ic = randi(numel(IC_mean.m));

    q0 = IC_mean.q(ind_ic,:);
    w0 = IC_mean.w(ind_ic,:);
    C0 = IC_mean.C(ind_ic,:);
    m0 = IC_mean.m(ind_ic,:);


end
