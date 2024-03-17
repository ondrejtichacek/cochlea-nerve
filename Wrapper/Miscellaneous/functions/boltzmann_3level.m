function [ y ] = boltzmann_3level( x, c1, c2, y1, y2 )
%BOLTZMANN_3LEVEL 

y = 1./(1 + c1*exp(-x./y1) + c2*exp(-x./y2));

end

