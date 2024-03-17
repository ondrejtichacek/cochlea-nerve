function dbprintf(s, varargin)
%DBPRINT
ST = dbstack(1);

% fprintf(['%s (%d): ', s], ST(1).name, ST(1).line, varargin{:});

cprintf(-[1,0.5,0], '%s (%d): ', ST(1).name, ST(1).line)
cprintf([1,0.5,0], s, varargin{:});


end

