function [fun, xbest, fval] = boltzmann(varargin)
              
[fun, xbest, fval] = fitting.optim_fit(@boltzmann, varargin{:});
        
end

