function [ m, h] = hh_mh_init( markovrates )
%HH_MH_INIT
%   Calculate initial number of open sodium channels according to Eqn. (7)
%   of Mino et al. (2002)

% Set initial values according to deterministic steady-state
m = markovrates.am ./ (markovrates.am + markovrates.bm);
h = markovrates.ah ./ (markovrates.ah + markovrates.bh);

end

