function [fun, xbest, fval] = gaussian(varargin)

[fun, xbest, fval] = fitting.optim_fit(@gaussian, varargin{:});
        
end

