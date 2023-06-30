function [theta] = Hill_Langmuir_A(L, n, KA)
%HILL_LANGMUIR 
% https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)
%
% theta ... fraction of the receptor bound by the ligand
% L     ... ligand concentration
% KA    ... ligand concentration producing half occupation
% n     ... Hill coeficient

theta = 1 ./ (1 + (KA./L).^n);

end

