function [ memLim ] = memoryLimit()
%MEMORYLIMIT 

st = dbstack;
if numel(st) < 2 || ~strcmp(st(2).file, 'memOpt.m')
    warning('Calling memoryLimit out of memOpt class.')
end

% Get free system memory
if ispc
    sts = memory();
    memLim = sts.MemAvailableAllArrays; % [Bytes]
elseif isunix
    sts = memoryunix();
    memLim = sts.AvailableMemory * 1000; % [Bytes]
elseif ismac
    error('not yet supproted');
end

end

