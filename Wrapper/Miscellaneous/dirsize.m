function [x] = dirsize(path)
s = dir(path);
name = {s.name};
isdir = [s.isdir] & ~strcmp(name,'.') & ~strcmp(name,'..');
subfolder = fullfile(path, name(isdir));
x = sum([s(~isdir).bytes cellfun(@dirsize, subfolder)]);
end