function [ parsedstats, allstats ] = memoryunix()
%MEMORYUNIX

if ~isunix
    error('Function MEMORYUNIX is not available on this platform.\n');
end

fid = fopen('/proc/meminfo', 'r');

mem = textscan(fid,'%s %f %s','delimiter',' ','multipleDelimsAsOne',true);    

fclose(fid);

names = mem{:,1};
values = mem{:,2};


names = strrep(names, ':','');
names = strrep(names, '(','_');
names = strrep(names, ')','');

for i = 1:length(names)
    allstats.(names{i}) = values(i);
end


parsedstats.AvailableMemory = allstats.MemFree + allstats.Buffers + allstats.Cached;

end

