function Phi = GHKflux_fast(V, Cin, Cout, P, z, T)
%GHKFLUX
% z ... charge [K]
% T ... temperature [K]
%
% Phi is the current density (flux) across the membrane carried by ion S, measured in amperes per square meter (A·m−2)
% arguments
%     V double % transmembrane potential in volts
%     Cin double % intracellular concentration of ion S, measured in mol·m−3 or mmol·l−1
%     Cout double % extracellular concentration of ion S, measured in mol·m−3 or mmol·l−1
%     P double = 1; % permeability of the membrane for ion S measured in m·s−1
% end

F = 96485.3328959;   % [C/mol]
R = 8.314459848;     % [J/mol/K]

% Phi = P * z^2 * (V*F^2 / (R*T)) * (Cin - Cout*exp(-z*V*F/(R*T))) / (1 - exp(-z*V*F/(R*T)));

U = z * V * F / R / T;
emU = exp(-U);
% epU = exp(U);

Phi = P .* z * F * U .* (Cin - Cout.*emU) ./ (1 - emU);

% Phi_in = P .* z * F * U .* Cin ./ (1 - emU);
% Phi_out = P .* z * F * U .* Cout ./ (1 - epU);

end

