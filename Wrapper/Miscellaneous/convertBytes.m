function [ res ] = convertBytes( varargin )
%CONVERTBYTES 

if nargin == 2
    from = to_Bytes(varargin{1});
elseif nargin == 3
    val = varargin{1};
    from = val*to_Bytes(sprintf('1 %s', varargin{2}));
end

to = to_Bytes(sprintf('1 %s', varargin{nargin}));

res = from/to;


end

