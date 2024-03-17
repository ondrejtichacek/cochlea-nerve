function [out, channels] = main(xgrid, mu, fs, args)
arguments
    xgrid (1,1) double
    mu (1,:) double
    fs (1,1) double
    
    args.V_hold = []
    args.I_hold = []
    
    args.EP = 90e-3 % Endocochlear potential
    args.EK = -70e-3 % K+ reversal potential
    args.Eextra = 0 %

    % args.Cm = 9.8e-12; %
    args.Cm = 12.0e-12; %
    
    args.channels = []

    args.random_MET (1,1) double = false
    args.num_MET_channels (1,1) double = 300
    
    args.random_K (1,1) double = false

    args.current_model (1,:) char = 'Ohm';
    % args.current_model (1,:) char = 'GHK';
    
    args.met_optim_params = []
end

% linear function f(x=0) = a, f(x=1) = b.
linf = @(a, b, x) a + (b-a) * x;

simulation_units = struct( ...
    'time', 's', ...
    'voltage', 'V', ...
    'current', 'A', ...
    'conductance', 'S', ...
    'capacitance', 'F', ...
    'charge', 'C');

%%

if isempty(args.channels)
    channels = BasolateralChannelsIHC('LopezPoveda_2006');
else
    channels = args.channels;
end

tauMet = 50e-6; % time constant of activation of MET channels

Gleak = 0e-9; % no membrane leakage
% Gleak = 2e-9;

% Gmet = 30e-9; %MET channel max conductance

% Cm = 9.8e-12;
% Cm = 12.0e-12;
Cm = args.Cm;

tauV = 1e-5; % time constant V_hold
tauCk = 10e-3; % time constant [K]
% tauCk = 0.1e-3; % time constant [K]

C_K_background = 4.8e-3;

% IkCK_fac = 1.34e6; % such that after equil, the steady state CK == Nernst_inv(EK, "conc_in", 131e-3);
%                   % = 3.1069e-3 for EK = -100
IkCK_fac = 12e6;

[~,~,~,~,MET_popen,Gmet_max] = IHC.METConductance(xgrid, simulation_units, ...
    optim_params=args.met_optim_params, ...
    version='Tichacek_2022');

popen = @(y) MET_popen(y * 1e9, ones(size(y))); % m to nm
Gmet = Gmet_max * 1e-9; % ns to S;

EP = args.EP;
EK = args.EK + args.Eextra;

% scale_by_tonotopy = true;
% if scale_by_tonotopy
%     GKf_x_fac = @(x) linf(1, 1, x);
%     GKs_x_fac = @(x) linf(1, 0.1, x);
%     GKn_x_fac = @(x) linf(0, 1, x);
%     
%     GKf = GKf * GKf_x_fac(xgrid);
%     GKs = GKs * GKs_x_fac(xgrid);
%     GKn = GKn * GKn_x_fac(xgrid);
% end

%%

dt = 1/fs;

% Numerical calculations based on the impulse invariance method
% https://en.wikipedia.org/wiki/Impulse_invariance

% coefficient for including activation time constants. 
alphaMet = exp(-dt./tauMet);
alphaCk = exp(-dt./tauCk);
alphaV = exp(-dt./tauV);

n_equil = round(100e-3*fs);
n_run = numel(mu);
n_GK = numel(channels);

out.Vm = zeros(n_run, 1);
out.Omet = zeros(n_run, 1);
out.Imet = zeros(n_run, 1);
out.Ileak = zeros(n_run, 1);
out.IK = zeros(n_run, n_GK);
out.CK = zeros(n_run, 1);
out.EK = zeros(n_run, 1);
out.GK = zeros(n_run, n_GK);

equil.Vm = zeros(n_equil, 1);
equil.Omet = zeros(n_equil, 1);
equil.Imet = zeros(n_equil, 1);
equil.Ileak = zeros(n_equil, 1);
equil.IK = zeros(n_equil, n_GK);
equil.CK = zeros(n_equil, 1);
equil.EK = zeros(n_equil, 1);
equil.GK = zeros(n_equil, n_GK);

if ~isempty(args.I_hold)

    if numel(args.I_hold) == 1
        args.I_hold = args.I_hold * ones(1, n_run);
    end

    args.I_hold = [args.I_hold(1) * ones(1,n_equil) args.I_hold];
end
if ~isempty(args.V_hold)

    if numel(args.V_hold) == 1
        args.V_hold = args.V_hold * ones(1, n_run);
    end

    args.V_hold = [args.V_hold(1) * ones(1,n_equil) args.V_hold];
end

mu = [zeros(1, n_equil) mu]; % input (cilia vibrations) padded with 100 ms of silence

mtIn = popen(mu);
mt = popen(0);

if args.random_MET == true
    for i = 1:length(mtIn)
        mtIn(i) = rand_open(mtIn(i), args.num_MET_channels);
    end
end


%%
% initial parameters
V = linf(-55.0e-3, -47.5e-3, xgrid); % approx resting potential

Ck = C_K_background;

order = [channels.order];
nC = numel(channels);
mK = zeros(1, nC);

%% 1st order channels

idx = order == 1;

channel_parameters = [channels(idx).parameters];

nC1 = numel(channel_parameters);

if nC1 > 0
    V0 = [channel_parameters.V];
    S0 = [channel_parameters.S];
    tau = [channel_parameters.tau];
    
    alphaK = exp(-dt./tau);
    
    GK_1 = [channel_parameters.G];
else
    GK_1 = [];
end
    

%% 2nd order channels

idx = order == 2;

channel_parameters = [channels(idx).parameters];

nC2 = numel(channel_parameters);

if nC2 > 0
    V1 = [channel_parameters.V1];
    V2 = [channel_parameters.V2];
    
    S1 = [channel_parameters.S1];
    S2 = [channel_parameters.S2];
    
    T1max = [channel_parameters.T1max];
    T2max = [channel_parameters.T2max];
    
    T1min = [channel_parameters.T1min];
    T2min = [channel_parameters.T2min];
    
    aT1 = [channel_parameters.aT1];
    aT2 = [channel_parameters.aT2];
    
    bT1 = [channel_parameters.bT1];
    bT2 = [channel_parameters.bT2];
    
    GK_2 = [channel_parameters.G];
else
    GK_2 = [];
end

%% Solve model

GK = [GK_1, GK_2];

if nC1 > 0
    mK(1:nC1) = 1./(1+exp(-(V-V0)./S0));
end

[MM, AA, ZZ] = deal([]);

j = nC1;
for i = 1:nC2
    j = j + 1;

    [~, ~, ZZ] = pOpenIHC_core(V, V1(i), S1(i), V2(i), S2(i), ...
        T1max(i), aT1(i), bT1(i), T1min(i), T2max(i), aT2(i), bT2(i), T2min(i), ...
        MM, AA, ZZ); % MM, AA and ZZ are passed just to speed up allocation

    P{i} = ZZ;
    mK(j) = P{i}(1);
end

GKt = 0*ones(1, numel(GK));

for it = 1:length(mtIn)

    % 1. MET current
    mt = alphaMet*mt + (1-alphaMet)*mtIn(it);
    Imet = Gmet*mt*(V - EP);

    % 2. Leakage current
    Ileak = Gleak*V;

    % 3. basolateral channel pOpen
    if nC1 > 0
        mK(1:nC1) = 1./(1+exp(-(V-V0)./S0)).*(1-alphaK) + mK.*alphaK;
    end

    j = nC1;
    for i = 1:nC2
        j = j + 1;

        [MM, AA, ZZ] = pOpenIHC_core(V, V1(i), S1(i), V2(i), S2(i), ...
            T1max(i), aT1(i), bT1(i), T1min(i), T2max(i), aT2(i), bT2(i), T2min(i), ...
            MM, AA, ZZ); % MM, AA and ZZ are passed just to speed up allocation

        dP = MM \ (-AA * P{i} + ZZ);

        P{i} = P{i} + dP*dt;
        
        mK(j) = P{i}(1);
    end

    % 4. Basolateral current
    switch args.current_model
        case 'Ohm'
            if args.random_K
                if it == 1
                    GKt = mK .* GK;
                end
                for l = 1:numel(GK)
                    GKt(l) = alphaMet*GKt(l) + (1-alphaMet) * GK(l) * rand_open(mK(l), round(GK(l)/15e-12));
                    % GKt(l) = GK(l) * rand_open(mK(l), round(GK(l)/15e-12));
                end

                Ik = GKt .* (V - EK);
            else
                Ik = mK .* GK .* (V - EK);
            end
        
        case 'GHK'
            C = Nernst_inv(EK, "conc_in", 131e-3);
            % C = Ck + C_K_background;
            Ik = potassium_current_GHK(V, GK, mK, C);
    end

    Ik_total = sum(Ik);

    % 5. Extracellular K concentration
    Ck = alphaCk*Ck + (1-alphaCk)*(Ik_total*IkCK_fac);

    % 6. Nernst potential
    EKN = Nernst_fast(131e-3, Ck + C_K_background, 1, 310);
    % EK = EKN;


    % 7. Total current
    I = Ileak + Imet + Ik_total;

    % 8. Current clamp
    if ~isempty(args.I_hold)
        I = I + args.I_hold(it);
        % I = args.I_hold(it);
    end

    % 9. Transmembrane voltage
    dV = -I/Cm;
    V = V + dV*dt;

    % 10. Voltage clamp
    if ~isempty(args.V_hold)
        % Vh = args.V_hold(it);
        % V = alphaV*V + (1-alphaV)*Vh;
        V = args.V_hold(it);
    end

    if it > n_equil
        jt = it - n_equil;
        out.Vm(jt) = V;
        out.Omet(jt) = mt;
        out.Ileak(jt) = Ileak;
        out.Imet(jt) = Imet;
        out.IK(jt,:) = Ik;
        out.CK(jt) = Ck + C_K_background;
        out.EK(jt) = EKN;
        % out.GK(jt,:) = mK;
        out.GK(jt,:) = GK .* mK;
    else
        jt = it;
        equil.Vm(jt) = V;
        equil.Omet(jt) = mt;
        equil.Ileak(jt) = Ileak;
        equil.Imet(jt) = Imet;
        equil.IK(jt,:) = Ik;
        equil.CK(jt) = Ck + C_K_background;
        equil.EK(jt) = EKN;
        equil.GK(jt,:) = GK .* mK;
    end
end

out.Itotal = out.Ileak + out.Imet + sum(out.IK, 2);


end

function conc_it()
  %% Calculate concentration
  
  ind_r_channel = 1;
  
  C_channels = zeros(numel(ps),1);
  if it > 1
    for i = 1:numel(ps)
      C_channels(i) = ps(i).concentration(ind_r_channel, it - 1);
    end
  end
  
  C_channels = C_channels + C_K_background;
  
  FAC = 1;

  [Ik, Ike, Ek] = potassium_current(V, GK/FAC, mK, C_channels);
  
  
  Ik_total = sum(Ik) * FAC;
  Ike_total = sum(Ike);

  for i = 1:numel(ps)
    ps(i).current(it) = Ike_total(i) / membrane_area; % minus because of convention
  end
  
  Ck = potassium_concentration(ps, it);

  for i = 1:numel(ps)
    ps(i).concentration(:, it) = Ck(:, i);
  end
end


function I = potassium_current_GHK(Vm, G, m, C)
arguments
  Vm (1,1) double
  G (:,1) double
  m (:,1) double % num_channel x 1
  C (:,1) double % num_channel x 1
end

% K+ ion charge
charge = 1;

% intracellular concentration
C_intracellular = 131e-3; % M

% calcium reversal potential
% E = Nernst(C_intracellular, C, 'charge', charge);

% amp_to_electron_per_second = 6.242e18;

% Ie = I * amp_to_electron_per_second / charge; % ions/sec

Cin = C_intracellular;
Cout = C;

area = 1;

me = 9.109383632e-31; % kg ... mass of electron
ce = -1.602e-19; % C ... charge of electron

fac = 3.9e-2*6.0e-15 ; % m^3 / s^2 ... whatever this is
P = G / ce^2 / charge^2 * me * fac;
P = P / area;

Cin = Cin*1000; % M to mol/m3
Cout = Cout*1000; % M to mol/m3
V = Vm;

P = m .* P;

[Phi, Phi_in, Phi_out] = GHKflux(V, Cin, Cout, P, charge = charge);

I = Phi * area;
% I_in = Phi_in * area;
% I_out = Phi_out * area;

end

function [C] = potassium_concentration(ps, it)
arguments
  ps
  it
end

N_A = 6.02214076e23; % mol^-1 (Avogadro constant)

C = zeros(ps(1).nr, numel(ps));
for i = 1:numel(ps)
  
  % j = ps(i).current(:, it);
  % current is handled internally by ps.iterate
    
  C(:,i) = ps(i).iterate(it);
end

C = C / N_A; % M

% c ... mol/m3
% convert c from m/m3 to mol/L = M (molar)
C = C / 1e3;

end
