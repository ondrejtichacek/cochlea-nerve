function [varargout] = findResultExtended(extra_dirs, PART, varargin)
arguments
    extra_dirs cell
    PART (1,:) char
end
arguments (Repeating)
    varargin
end
%FINDRESULTEXTENDED 

n = max([1, nargout]);

[varargout{1:n}] = findResult(PART, varargin{:});

for i = 1:numel(extra_dirs)
    runopt = varargout{1};
    
    if runopt.found.(PART) == false
        
        [varargout{1:n}] = findResult(PART, varargin{:}, ...
            'result_base_dir', extra_dirs{i});
    else
        break
    end
end

end

