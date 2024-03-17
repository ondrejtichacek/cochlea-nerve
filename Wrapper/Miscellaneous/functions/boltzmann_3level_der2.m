function [ ddy ] = boltzmann_3level_der2( x, c1, c2, y1, y2 )
%BOLTZMANN_3LEVEL 

%y = 1./(1 + c1*exp(-x./y1) + c2*exp(-x./y2));

% second derivative
% (2 (-(c1 e^(-x/y1))/y1 - (c2 e^(-x/y2))/y2)^2)/(c1 e^(-x/y1) + c2 e^(-x/y2) + 1)^3 - ((c1 e^(-x/y1))/y1^2 + (c2 e^(-x/y2))/y2^2)/(c1 e^(-x/y1) + c2 e^(-x/y2) + 1)^2
ddy = (2*(-(c1 *exp(-x/y1))/y1 - (c2 *exp(-x/y2))/y2)^2)/(c1 *exp(-x/y1) + c2 *exp(-x/y2) + 1)^3 - ((c1 *exp(-x/y1))/y1^2 + (c2 *exp(-x/y2))/y2^2)/(c1 *exp(-x/y1) + c2 *exp(-x/y2) + 1)^2;

end

