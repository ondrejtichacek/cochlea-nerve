function [MM, AA, ZZ] = pOpenIHC_core( ...
    V, V1, S1, V2, S2, T1max, aT1, bT1, T1min, T2max, aT2, bT2, T2min, MM, AA, ZZ)
%POPENIHC_CORE

% Numstacks = numel(V);
% n = 2*Numstacks;
% 
% assert(Numstacks == 1)

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

% calculate voltage dependent time constatns tau1 and tau2
tau1 = T1min + (T1max - T1min)./(1 + exp((aT1 + V)/bT1));
tau2 = T2min + (T2max - T2min)./(1 + exp((aT2 + V)/bT2));

% and the voltage dependent steady-state open probability
Oinf = 1./(1 + exp((V1 - V)/S1) .* (1 + exp((V2 - V)/S2)));

% assemble the RHS
% assembling the full system is somewhat trickier:

if ~isempty(MM) % for speed
    ZZ(2) = Oinf;
    AA(2,2) = tau1+tau2;
    MM(2,2) = tau1.*tau2;
else
    ZZ = [0; Oinf];
    MM = [1 0; 0 tau1.*tau2];
    AA = [0 -1; 1 tau1+tau2];    
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

%% Currently unused

% dtau1 = (T1min - T1max)./(2*bT1*cosh((aT1 + V)/bT1) + 2*bT1);
% dtau2 = (T2min - T2max)./(2*bT2*cosh((aT2 + V)/bT2) + 2*bT2);
% 
% dP_inf = ( ...
%         exp((V1*S2 + V*(S2 + S1))/(S2*S1)) ...
%         .* ...
%         (exp(V2/S2)*(S2 + S1) + S2*exp(V/S2)) ...
%          ) ...
%     ./(S2*S1*( ...
%         exp(V1/S1 + V2/S2) ...
%       + exp(V1/S1 + V/S2) ...
%       + exp(V*(1/S2 + 1/S1)) ...
%              ).^2);
% 
% JJ = -AA;
% 
% dP_V = -(dtau1 + dtau2);
% dZ_V = dP_inf;

end