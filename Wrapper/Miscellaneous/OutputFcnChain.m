function status = OutputFcnChain(fun, varargin)
status = zeros(1, numel(fun));
for i = 1:numel(fun)
    status(i) = fun{i}(varargin{:});        
end
status = any(status);
end