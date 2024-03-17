function [fun, xbest, fval] = boltzmann_der(varargin)

[fun, xbest, fval] = fitting.optim_fit(@boltzmann_der, varargin{:});

end