function [ varargout ] = disp_toc( t, SPEC, UNIT )
%DISP_TOC 

if nargin < 1
    t = toc;
else
    assert(isnumeric(t), 'The first argument must be a number')
    if isa(t, 'uint64')
        warning('The class of the first argument is uint64. See if you are calling by mistake an equivalent of t = tic, disp_toc(t).')
    end
end

if nargin < 2
    
    if abs(t) < 0.001
        
        STR = sprintf('%g %s', t*1e6, 'usec');

    elseif abs(t) < 1
        
        STR = sprintf('%.3g %s', t*1e3, 'msec');                
        
    elseif abs(t) < 60
        
        STR = sprintf('%.3g %s', t, 'sec');
        
    elseif abs(t) < 3600
        
        SPEC = 'MM:SS';
        UNIT = 'minutes';
        
        STR = sprintf('%s %s', datestr(t/86400, SPEC), UNIT);
        
    elseif abs(t) < 86400
        
        SPEC = 'HH:MM:SS';
        UNIT = 'hours';
        
        STR = sprintf('%s %s', datestr(t/86400, SPEC), UNIT);
        
    else
        
        d = floor(t/86400);
        h = mod(t, 86400);
        
        if d == 1 
            s = 'day';
        else
            s = 'days';
        end 
        
        SPEC = 'HH:MM:SS';
        UNIT = 'hours';
        
        STR = sprintf('%d %s %s %s', d, s, datestr(h/86400, SPEC), UNIT);
    end
    
else
    
    STR = sprintf('%s%s', datestr(t/86400, SPEC), UNIT);
    
end

if nargout == 0
    fprintf('Elapsed time is %s.\n',STR);
    
elseif nargout == 1
    varargout = {STR};
    
elseif nargout == 2
    varargout = {STR, t};
    
else
    error('Too many output arguments');
    
end

end
