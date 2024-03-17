function deleteResultFiles( FOLDER )
%DELETERESULTFILES 

d =  dir(fullfile(FOLDER));

names = {d.name};
is_dir = [d.isdir];

[names, I] = setdiff(names, {'.','..'});        
is_dir = is_dir(I);

for i = 1:length(names)
    name = names{i};

    if is_dir(i)
        srmdir(fullfile(FOLDER, name), 's');
    else
        sdelete(fullfile(FOLDER, name));
    end                
end

if ~isempty(names)
    warning('Result files removed:\n%s\n%s', ...
        FOLDER, ...
        sprintf('    %s\n', names{:}));
end

end
