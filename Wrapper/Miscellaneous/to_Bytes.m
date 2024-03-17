function val = to_Bytes(str)
% Parses string to number of bytes
%
% example:
%   to_Bytes('10 kB')  = 10000
%   to_Bytes('10 kiB') = 10240
%   to_Bytes('10 kb')  = 1250
%   to_Bytes('10 kib') = 1280
%   to_Bytes('1 MB')   = 1000000
%   to_Bytes('1 MiB')  = 1048576


C = strsplit(strtrim(str));

assert(numel(C) == 2);

numstr = C{1};
unitstr = C{2};

power = 0;
base = 1000;

len = length(unitstr);
assert(len <= 3, 'length(unit) must be <= 3; unit = `%s`.', unitstr)

unit = bit_or_byte(unitstr(len));

if len >= 2
    power = unit_power(unitstr(1));
end
if len == 3
    if unitstr(2) == 'i'
        base = 1024;
    else
        error('Undefined base `%s`.', unitstr(2))
    end
end

val = str2double(numstr) * base^power*unit;

    function val = bit_or_byte(b)
        % return multiplicator from `b` to bytes

        if b == 'b'
            val = 1/8; % bit
        elseif b == 'B'
            val = 1; % byte
        end
    end
    function val = unit_power(p)
        % return power of unit with the base of 1000 or 1024

        switch upper(p)
            case 'K', val = 1;
            case 'M', val = 2;
            case 'G', val = 3;
            case 'T', val = 4;
            case 'P', val = 5;
            case 'E', val = 6;
            case 'Z', val = 7;
            case 'Y', val = 8;
            otherwise
                error('Unknown power `%s`.', p)
        end
    end
end