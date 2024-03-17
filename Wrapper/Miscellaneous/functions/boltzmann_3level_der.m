function [ dy ] = boltzmann_3level_der( x, c1, c2, y1, y2 )
%BOLTZMANN_3LEVEL_DER

%y = 1./(1 + c1*exp(-x./y1) + c2*exp(-x./y2));

% derivative
% -(-(c1 e^(-x/y1))/y1 - (c2 e^(-x/y2))/y2)/(c1 e^(-x/y1) + c2 e^(-x/y2) + 1)^2
% alt: (c1 y2 e^(x (1/y1 + 2/y2)) + c2 y1 e^(x (2/y1 + 1/y2)))/(y1 y2 (c1 e^(x/y2) + c2 e^(x/y1) + e^(x (1/y1 + 1/y2)))^2)
% alt: (e^(x/y1 + x/y2) (c1 y2 e^(x/y2) + c2 y1 e^(x/y1)))/(y1 y2 (c1 e^(x/y2) + c2 e^(x/y1) + e^(x/y1 + x/y2))^2)

dy = -(-(c1 *exp(-x/y1))/y1 - (c2 *exp(-x/y2))/y2)./(c1 *exp(-x/y1) + c2 *exp(-x/y2) + 1).^2;

end

