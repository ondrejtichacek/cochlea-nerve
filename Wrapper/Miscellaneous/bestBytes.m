function [ val, unit ] = bestBytes( value, unit_in, base )
%BESTBYTES
% Returns 'value' of 'unit_in' converted to 'unit' of base 'base' such that
%   1.   value / base^k < base
%   2.   value / base^k >= 1
%
%   => k = floor(log(value)/log(base))
%
%   3.   k is truncuated so that k >= 0 and k <= K, where K = 8
%

if nargin == 1
    if ischar(value)
        value = to_Bytes(value);
    end
elseif nargin == 2
    value = convertBytes(value, unit_in, 'B');
end

if nargin < 3
    base = 1000;
end

k = floor(log(value)/log(base));

if k <= 0
    unit = 'B';
else
    
    if base == 1024
        b = 'i';
    else
        b = '';
    end
    
    U = {'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};
    
    unit = [ U{min([k,numel(U)])}, b, 'B' ];
end

val = convertBytes(value, 'B', unit);

end
