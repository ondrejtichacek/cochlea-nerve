function V = Nernst_fast(Cin, Cout, z, T)
%NERNST
% z ... charge [K]
% T ... temperature [K]

F = 96485.3328959;   % [C/mol]
R = 8.314459848;     % [J/mol/K]


V = R.*T./z./F.*log(Cout./Cin); % [J/mol/K*K/C*mol] = [J/C] = [V]

end