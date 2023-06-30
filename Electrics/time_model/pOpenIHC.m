function [MM, AA, ZZ, JJ, dP_V, dZ_V, MM_V_pat, Oinf] = pOpenIHC(V, channel, do)
%POPENIHC 
arguments
    V (:,1) double
    channel (1,1) struct
    do.mass (1,1) logical = true
    do.system (1,1) logical = true
    do.rhs (1,1) logical = true
    do.jacobian (1,1) logical = true
end

Numstacks = numel(V);
n = 2*Numstacks;

persistent spdiag_N
if isempty(spdiag_N)
    spdiag_N = spdiags(ones(Numstacks,1),0,Numstacks,Numstacks);
end

% G = channel.parameters.G;          % Maximum conductance
% we don't use the conductance here, it is used in the circuit 

V1 = channel.parameters.V1;        % Half-activation setpoint
S1 = channel.parameters.S1;        % Voltage sensitivity constant
V2 = channel.parameters.V2;        % Half-activation setpoint
S2 = channel.parameters.S2;        % Voltage sensitivity constant

% Time constants
T1max = channel.parameters.T1max;
aT1 = channel.parameters.aT1;
bT1 = channel.parameters.bT1;
T1min = channel.parameters.T1min;

T2max = channel.parameters.T2max;
aT2 = channel.parameters.aT2;
bT2 = channel.parameters.bT2;
T2min = channel.parameters.T2min;

%%
% We have a 2nd order ODE
% t1 t2 O'' + (t1 + t2) O' + O - Oinf = 0
%
% where
%        O ... channel open probability
%    t1,t2 ... time constants
%     Oinf ... steady state open probability
%
% Because its a second order ODE, we need to transform it to 2 first order
% ODEs as follows:
% 
%   P2 = O'
%   P1 = O
%
%   where P1' = P2
%
% And we have
%
%         P1' - P2 = 0
%   t1 t2 P2' + (t1+t2) P2 + P1 - Oinf = 0
%
% or equiv.
%
%         P1' = P2
%   t1 t2 P2' = -((t1+t2) P2 + P1) + Oinf = 0
%
% We can write this in a matrix form as
%
% M*dP/dt = -A*P + Z =: f(t)
%
% M = [ 1     0
%       0   tau1.*tau2 ];
% 
% A = [ 0     -1
%       1  tau1+tau2 ];
%
% Z = [ 0 
%       Oinf ];
%
% The sought solution open probability O is therefore sol. of eq. 1
%

%% System

if do.mass || do.system || do.jacobian
    % calculate voltage dependent time constatns tau1 and tau2
    tau1 = T1min + (T1max - T1min)./(1 + exp((aT1 + V)/bT1));
    tau2 = T2min + (T2max - T2min)./(1 + exp((aT2 + V)/bT2));

    if false
        VV = 1e-3*(-100:0.1:0);
        pl_tau1 = T1min + (T1max - T1min)./(1 + exp((aT1 + VV)/bT1));
        pl_tau2 = T2min + (T2max - T2min)./(1 + exp((aT2 + VV)/bT2));
        figure
        hold on
        plot(VV*1e3, pl_tau1*1e3, 'DisplayName', '\tau_1')
        plot(VV*1e3, pl_tau2*1e3, 'DisplayName', '\tau_2')
        legend()
        xlabel('voltage (mV)')
        ylabel('\tau (ms)')
    end
end
if do.rhs
    % and the voltage dependent steady-state open probability
    Oinf = 1./(1 + exp((V1 - V)/S1) .* (1 + exp((V2 - V)/S2)));

    assert(size(Oinf,2) == 1)
    assert(size(Oinf,1) == Numstacks)

    % assemble the RHS
    ZZ = [zeros(1, Numstacks); Oinf(:)'];
    % this a bit of hack
    ZZ = ZZ(:);
    assert(ZZ(1) == 0)
else
    ZZ = [];
    Oinf = [];
end

% assembling the full system is somewhat trickier:
if do.mass
    if Numstacks == 1
        MM = [1 0; 0 tau1.*tau2];
    else
        % the mass could be generally assembled like this
        
        M11 = [1 0; 0 0];
        M22 = [0 0; 0 1];
    
        MM11 = kron(spdiag_N, M11);
        % MM11 = kron(spdiags(ones(Numstacks,1),0,Numstacks,Numstacks), M11);
        % MM22 = kron(spdiag_N .*(tau1.*tau2), M22); % slower here
        MM22 = kron(spdiags(tau1.*tau2,0,Numstacks,Numstacks), M22);
    
        MM = MM11 + MM22;
    
        MM_V_pat = (MM22 ~= 0); % voltage dependent part of MM;
        
        % or like this (actually slower)
        
        % mm = [ones(1,Numstacks); tau1'.*tau2'];
        % MM = spdiags(mm(:),0,n,n);
        
        % mm = [zeros(1,Numstacks); ones(1,Numstacks)];
        % MM_V_pat = spdiags(mm(:),0,n,n); % voltage dependent part of MM;
    end
else
    [MM, MM_V_pat] = deal([]);
end

if do.system || do.jacobian
    if Numstacks == 1
        AA = [0 -1; 1 tau1+tau2];
    else
        AA = spalloc(n, n, 3*Numstacks); % sparse zero
    
        nzmaxAA = nzmax(AA);
    
        A1221 = [0 -1; 1 0];
        A22 =   [0  0; 0 1];
        AA = AA + kron(spdiag_N, A1221);
        % AA = AA + kron(spdiags(ones(Numstacks,1),0,Numstacks,Numstacks), A1221);
        % AA = AA + kron(spdiag_N .* (tau1+tau2), A22); % slower here
        AA = AA + kron(spdiags(tau1+tau2,0,Numstacks,Numstacks), A22);
    
        assert(nzmax(AA) == nzmaxAA)
    end
else
    AA = [];
end
    
% MM*dPP/dt = -AA*PP + ZZ

%% Jacobian
% Calculate Jacobian
%
%  J = df/dP
%
% via the product rule
%  J = -A(t,P) - dA(t,P)/dP * P + dZ(t,P)/dP
%
% What we actually have
%  J = -A(t,P) - dA(t,V)/dY * Y + dZ(t,V)/dY
%
%
if do.jacobian
    dtau1 = (T1min - T1max)./(2*bT1*cosh((aT1 + V)/bT1) + 2*bT1);
    dtau2 = (T2min - T2max)./(2*bT2*cosh((aT2 + V)/bT2) + 2*bT2);

    if false
        VV = 1e-3*(-100:0.1:0);
        pl_dtau1 = (T1min - T1max)./(2*bT1*cosh((aT1 + VV)/bT1) + 2*bT1);
        pl_dtau2 = (T2min - T2max)./(2*bT2*cosh((aT2 + VV)/bT2) + 2*bT2);
        figure
        hold on
        plot(VV, pl_dtau1, 'DisplayName', 'd\tau_1/dV')
        plot(VV, pl_dtau2, 'DisplayName', 'd\tau_2/dV')
    end
    dP_inf = ( ...
            exp((V1*S2 + V*(S2 + S1))/(S2*S1)) ...
            .* ...
            (exp(V2/S2)*(S2 + S1) + S2*exp(V/S2)) ...
             ) ...
        ./(S2*S1*( ...
            exp(V1/S1 + V2/S2) ...
          + exp(V1/S1 + V/S2) ...
          + exp(V*(1/S2 + 1/S1)) ...
                 ).^2);

    if any(isnan(dP_inf))
        ind = isnan(dP_inf);
        if V(ind) > 100e-3 % 100mV
            % fix limit for V*S -> infty
            disp(V(ind))
            dP_inf(ind) = 0;
        else
            error("not fixed")
        end
    end

    dtau1 = dtau1(:);
    dtau2 = dtau2(:);

    dP_inf = dP_inf(:)';

    assert(size(dP_inf,1) == 1)
    assert(size(dP_inf,2) == Numstacks)

    JJ = -AA;

    dP_V = -(dtau1 + dtau2);
    dZ_V = dP_inf;
else
    [JJ, dP_V, dZ_V] = deal([]);
end
end