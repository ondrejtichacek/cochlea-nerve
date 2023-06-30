function [ dz ] = TransductionRHS_v5( t, z, ...
    size_info, channels, vesicles, ps, rate_y, rate_l, rate_x, rate_r, G_Ca, C_Ca_background, ...
    channels_open_ss_parameters_normal, channels_open_ss_parameters_burst, ...
    transmitter_release_parameters, ...
    V_steady_state, ...
    dt, Vt)
arguments
    t
    z
    size_info
    channels
    vesicles
    ps
    rate_y
    rate_l
    rate_x
    rate_r
    G_Ca
    C_Ca_background
    channels_open_ss_parameters_normal
    channels_open_ss_parameters_burst
    transmitter_release_parameters
    V_steady_state
    dt
    Vt
end
%TRANSDUCTIONRHS

Vt = Vt(:);

n = size(Vt,1);

m = decompose_z(z, 'CaV13_channels_fraction', size_info);
Ca_blocked = decompose_z(z, 'CaV13_channels_blocked', size_info);
C_vesicles = decompose_z(z, 'Ca_concentration', size_info);
I = decompose_z(z, 'Ca_current', size_info);
q = decompose_z(z, 'NT_free', size_info);
c = decompose_z(z, 'NT_cleft', size_info);
w = decompose_z(z, 'NT_reprocessing', size_info);

c_proton = decompose_z(z, 'proton_cleft', size_info);

m_old = m;
I_old = I;
Ca_blocked_old = Ca_blocked;
C_old = C_vesicles;

% iteration index
it = ceil(t/dt);

%%

[d_state_inactivated, d_state_normal, d_state_burst] = deal(zeros(1,n));

% number of Ca_V1.3 channels in the vicinity of the synapse
num_CaV13 = channels.num;


%% mode switching
% Currently not used, see TransductionRHS_v4

%%

p_block = CaV13.CaProtonBlock(c_proton);

% assert(dt / channels.tau <= 1);
assert(dt / channels.tau_blocked <= 1);

%% Calculate rates

% S0t = 5.9e-3;
% V0t = -20.6e-3;
% S0t = 7e-3;
% V0t = -30e-3;
S0t = channels.S0t;
V0t = channels.V0t;

% alpha = -60; % v1
% alpha = -130; % v2
alpha = channels.alpha;
beta = alpha + 1/S0t;

% beta = -20;
% alpha = beta - 1/S0t;

% kp0 = 1/(3e-3); % v1
% kp0 = 1/(1.5e-5); % v2
kp0 = channels.kp0;

km0 = kp0 * exp(V0t * (beta - alpha));

% S0 = 1 / (beta - alpha); % V
% V0 = log(km0/kp0) / (beta - alpha); % V


%% compute conductance

kp = kp0 * exp(-alpha * Vt);
km = km0 * exp(-beta * Vt);

Q = [1 - kp*dt, km *dt; 
     kp*dt, 1 - km*dt];

cQ = cumsum(Q,1);

ST = ['c', 'o'];

for ii = 1:channels.num

    if channels.state(ii) == 'i' % inactivated
        continue
    end

    % unblock channel
    if channels.state(ii) == 'b' % blocked
        if t >= channels.tblocked(ii)
            channels.state(ii) = 'c';
        else
            continue
        end
    end

    % block channel
    if rand(1) < p_block * dt / channels.tau_blocked

        % generate random close time from a gamma distribution
        % such that the mean (= a*b) is equal to the channels.tau_blocked
        %     a ... shape
        %     b ... scale

        a = 4;
        b = channels.tau_blocked / a;
        % tau = gamrnd(a, b, 1, 1);
        tau = b * gammaincinv(rand(1,1), a); % gamma random number without toolbox

        channels.state(ii) = 'b';
        channels.tblocked(ii) = t + tau;

    else
        r = rand(1);
    
        if channels.state(ii) == 'c' % closed
            st = 1;
        elseif channels.state(ii) == 'o' % open
            channels.topen(ii) = t; % last open time
            st = 2;
        else
            error('unknown state %s', channels.state(ii))
        end
        
        st_new = find(r < cQ(:,st), 1);
    
        channels.state(ii) = ST(st_new);

    end
end


m = sum(channels.state == 'o') / num_CaV13;
Ca_blocked = sum(channels.state == 'b') / num_CaV13;

%% Calculate concentration
% currently we calculate at two points:
%   1. close to the vesicle
%   2. close to the ion channel

all_channel_switch = true;
log_partial_concentrations = false;

if all_channel_switch
    ind_r_channel = 1:numel(ps);
else
    ind_r_channel = 1;
end

C_channels = zeros(numel(ps),1);
if it > 1
    for i = 1:numel(ps)
        if all_channel_switch
            if log_partial_concentrations
                C = ps(i).concentration(ind_r_channel, it - 1);
            else
                C = ps(i).concentration(ind_r_channel);
            end
            C_channels(i) = sum(C);
            % ratio(i) = 1 - C(i) / sum(C);
        else
            C_channels(i) = ps(i).concentration(ind_r_channel, it - 1);
        end
    end
    % disp(100*mean(ratio));
end

C_channels = C_channels + C_Ca_background;

open_channels = channels.state == 'o';

% I_ohm = calcium_current(Vt, G_Ca, open_channels, C_channels);
I_GHK = calcium_current_GHK(Vt, G_Ca, open_channels, C_channels);

% disp(mean(I_GHK(open_channels) ./ I_ohm(open_channels)));

I = I_GHK;

for i = 1:numel(ps)
    ps(i).current(it) = - I(i); % minus because of convention
    ps(i).lastopen = channels.topen(i);
end

[C_vesicles, C_all] = calcium_concentration(vesicles, ps, it, all_channel_switch);


for i = 1:numel(ps)
    if log_partial_concentrations
        ps(i).concentration(:, it) = C_all(:, i);
    else
        ps(i).concentration = C_all(:, i);
    end
end

C_vesicles = C_vesicles + C_Ca_background;

% vesicles.Ca_concentration(:, it) = C_vesicles;

%% Transmitter Release and Recycling

[dq, dc, dw, dc_proton, vesicles] = NTdynamicsRHS_v5_core( t, ...
    q, c, w, c_proton, ...
    vesicles, rate_y, rate_l, rate_x, rate_r, ...
    transmitter_release_parameters, ...
    dt, C_vesicles);

%% Build dz

dm = (m - m_old) / dt;

dCa_blocked = (Ca_blocked - Ca_blocked_old) / dt;

charge = 2;
amp_to_electron_per_second = 6.242e18;

I = sum(I) / amp_to_electron_per_second * charge;
dI = (I - I_old) / dt;

dC = (C_vesicles - C_old) / dt;

dz = [dm; dCa_blocked; dI; dC; dq; dc; dw; dc_proton];

d_state = [d_state_inactivated; d_state_normal; d_state_burst] / dt;

dz = [dz; d_state];

if any(isnan(dz))
    error('NaN encountered')
end

end

function I = calcium_current(Vm, G, m, C)
arguments
    Vm (1,1) double
    G (1,1) double
    m (:,1) double  % num_channel x 1
    C (:,1) double  % num_channel x 1
end

% Ca2+ ion charge
charge = 2;

% extracellular concentration
C_extracellular = 1.3e-3; % M

% calcium reversal potential
E = Nernst(C, C_extracellular, 'charge', charge); % this is something new !!!

amp_to_electron_per_second = 6.242e18;

I = m .* G .* (Vm - E); % A

I = I * amp_to_electron_per_second / charge; % ions/sec

% figure
% plot(t*1e3, chi.*j/1e3);
% ylabel('Ca2+ flux (ion/ms)')
% xlabel('time (ms)')

end

function I = calcium_current_GHK(Vm, G, m, C)
arguments
    Vm (1,1) double
    G (1,1) double
    m (:,1) double  % num_channel x 1
    C (:,1) double  % num_channel x 1
end

% Ca2+ ion charge
charge = 2;

% extracellular concentration
C_extracellular = 1.3e-3; % M

T = 310; % [K] ~ 37C

area = 1;

me = 9.109383632e-31; % kg ... mass of electron
ce = -1.602e-19; % C ... charge of electron

fac = 6.0e-15 ; % m^3 / s^2 ... whatever this is
P = G .* m / ce^2 / charge^2 * me * fac;
P = P / area;

Cin = C*1000; % M to mol/m3
Cout = C_extracellular*1000; % M to mol/m3

Phi = GHKflux_fast(Vm, Cin, Cout, P, charge, T);

I = Phi * area;

amp_to_electron_per_second = 6.242e18;

I = I * amp_to_electron_per_second / charge; % ions/sec

% figure
% plot(t*1e3, chi.*j/1e3);
% ylabel('Ca2+ flux (ion/ms)')
% xlabel('time (ms)')

end

function [C_vesicle, C] = calcium_concentration(vesicles, ps, it, all_channel_switch)
arguments
    vesicles
    ps
    it
    all_channel_switch
end

N_A = 6.02214076e23; % mol^-1 (Avogadro constant)

C = zeros(ps(1).nr, numel(ps));
for i = 1:numel(ps)
    
    % j = ps(i).current(:, it);
    % current is handled internally by ps.iterate
        
    C(:,i) = ps(i).iterate(it);
end


if all_channel_switch
    ind_r_channel = numel(ps);
else
    ind_r_channel = 1;
end

C_vesicle = zeros(vesicles.num, 1);
for jj = 1:vesicles.num
    
    % iterate over channels in vicinity of the vesicle
    close_channels = vesicles.close_channels{jj};
    C_vesicle(jj) = sum(C(ind_r_channel + jj, close_channels));
    
end

C = C / N_A; % M
C_vesicle = C_vesicle / N_A; % M

% c ... mol/m3
% convert c from m/m3 to mol/L = M (molar)
C = C / 1e3;
C_vesicle = C_vesicle / 1e3;

end

function v = decompose_z(z, variable, size_info)

si = size_info.(variable);
v = z(si.start : si.end);

end
