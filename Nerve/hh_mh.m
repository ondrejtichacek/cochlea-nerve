function [ dm, dh ] = hh_mh( m, h, markovrates )
%HH_MH
%   Calculate number of open sodium channels according to Eqns. (7) & (8)
%   of Mino et al. (2002)

dm = (markovrates.am .* (1-m)) - (markovrates.bm .* m);
dh = (markovrates.ah .* (1-h)) - (markovrates.bh .* h);

end

