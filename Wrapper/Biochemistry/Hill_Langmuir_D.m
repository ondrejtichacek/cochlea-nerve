function [theta] = Hill_Langmuir_D(L, n, KD)
%HILL_LANGMUIR 
% https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)
%
% theta ... fraction of the receptor bound by the ligand
% L     ... ligand concentration
% KD    ... apparent dissociation constant derived from the law of mass action
% n     ... Hill coeficient

% theta = L.^n ./ (KD + L.^n);
theta = 1 ./ (1 + KD./(L.^n));

end

