function [ dy ] = derHH( t, y, I_input, t_input, Na_max, Rm, gNa, E_Na, Cm, calcN_Na, mh, dt, I )
%DERHH 

% Stimulating current
if nargin < 13
    I = interp1( t_input, I_input, t );
end

I = I(:);

I(isnan(I)) = 0;

n = size(I,1);
    
% Extract individual variables
V = y(1:n);
m = y(n+(1:n));
h = y(2*n+(1:n));

% Calculate sodium gating particle transition rates for present
% membrane potential
markovrates = calcrates(V);

[dm, dh] = mh(m, h, markovrates);

N_Na = calcN_Na(m, h, markovrates, Na_max);

% Calculate membrane potential for next Euler time step
dV = (I - V/Rm - gNa * N_Na .* (V - E_Na)) / Cm;

dy = [dV; dm; dh];
% dy = Frequency([dV; dm; dh], 'kHz');

end

