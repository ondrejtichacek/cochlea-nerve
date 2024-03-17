function N_Na = inithhfloor(markovrates, Na_max, dt)
%INITHHFLOOR
%   Calculate initial number of open sodium channels according to Eqn. (7)
%   of Mino et al. (2002) with floor(N_Na)

global m h;
    
% Set initial values according to deterministic steady-state
m = markovrates.am ./ (markovrates.am + markovrates.bm);
h = markovrates.ah ./ (markovrates.ah + markovrates.bh);

N_Na = floor(Na_max * m.^3 .* h);
    
end

