function N_Na = hhfloor(markovrates, Na_max, dt)
%HHFLOOR
%   Calculate number of open sodium channels according to Eqns. (7) & (8)
%   of Mino et al. (2002) with floor(N_Na)

global m h;

m = m + ((markovrates.am .* (1-m))-(markovrates.bm .* m))*dt;
h = h + ((markovrates.ah .* (1-h))-(markovrates.bh .* h))*dt;

N_Na = floor(Na_max * m.^3 .* h);

end

