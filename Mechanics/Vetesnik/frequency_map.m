function [varargout] = frequency_map(args)
%FREQUENCY_MAP
arguments
    args.frequency double = []
    args.position double = []
end

assert(~(~isempty(args.frequency) && ~isempty(args.position)), ...
    'only one of `frequency` and `position` can be specified')

%%

base = 2.4531;
alpha = -6.36;
Fo = 21100;

char_frequency = @(x) Fo * base.^(alpha*(x));
char_position = @(c) log(c/Fo)/log(base^alpha);

%%

if ~isempty(args.frequency)
    
    x = char_position(args.frequency);
    
    if any(x > 1)
        error('position overflow (x > 1)')
    elseif any(x < 0)
        error('position underflow (x > 1)')
    end
    
    varargout = {x};
    
elseif ~isempty(args.position)
    
    x = args.position;
    
    if any(x > 1)
        error('position overflow (x > 1)')
    elseif any(x < 0)
        error('position underflow (x > 1)')
    end
    
    f = char_frequency(args.position);
    
    varargout = {f};
    
else
    varargout = {char_frequency, char_position};
    
end

end

