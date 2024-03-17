function [ S ] = matfileWhosByName( FILE, NAME, STAT, OPT )
%MATFILEWHOSBYNAME 

if nargin < 4
    OPT = {};
end
if nargin < 3
    STAT = [];
end
if nargin < 2
    NAME = [];
end


if nargin >= 1
    switch class(FILE)
        case 'char'
            W = whos('-file', FILE);
        case 'matlab.io.MatFile'
            W = whos(FILE);
        otherwise
            error('First argument must be either file name or matfile handle');
    end
    
    if isempty(W)
        S = [];
        return
    end
end

if isempty(NAME) && isempty(STAT)
    for i = 1:length(W)
        S.(W(i).name) = W(i);
    end
    
elseif isempty(NAME) && ~isempty(STAT)
    
    for i = 1:length(W)
        S.(W(i).name) = W(i).(STAT);
    end
    
elseif ~isempty(NAME) && isempty(STAT)
    
    [~,i] = find(strcmp({W.name},NAME),1);
    
    if ~isempty(i)
        S = W(i);
    else
        error('Variable %s not found in the file', NAME)
    end
    
else
    
    [~,i] = find(strcmp({W.name},NAME),1);
    if ~isempty(i)
        S = W(i).(STAT);
    else
        if any(strcmp('notFoundNoError', OPT))
            S = [];
        else
            error('Variable %s not found in the file', NAME)
        end
    end
    
end


end
