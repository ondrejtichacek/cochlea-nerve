function markovrates = calcrates(V)
%CALCRATES
%   Calculate sodium gating particle transition rates according to Eqns.
%   (2)-(5) of Mino et al. (2002)

markovrates.am = (1.872*(V - 25.41)) ./ (1 - exp((25.41 - V)/6.06));
markovrates.bm = (3.973*(21.001 - V)) ./ (1 - exp((V - 21.001)/9.41));
markovrates.ah = (-0.549*(27.74 + V)) ./ (1 - exp((V + 27.74)/9.06));
markovrates.bh = (22.57) ./ (1 + exp((56.0 - V)/12.5));

end

