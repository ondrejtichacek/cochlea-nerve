function N_Na = hhnint( m, h, markovrates, Na_max )
%HHNINT
%   Calculate number of open sodium channels according to Eqns. (7) & (8)
%   of Mino et al. (2002) with nint(N_Na)

N_Na = round(Na_max * m.^3 .* h);

end