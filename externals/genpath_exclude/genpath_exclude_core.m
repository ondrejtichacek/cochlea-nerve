function p = genpath_exclude_core(d, exclude_dir_str, exclude_dir_with_file_str)
arguments
    d (1,:) char
    exclude_dir_str char = ''
    exclude_dir_with_file_str char = ''
end
    p = '';

    % Generate path based on given root directory
    contents = dir(d);

    % set logical vector for subdirectory entries in d
    isdir = logical(cat(1,contents.isdir));
    
    if ~isempty(exclude_dir_with_file_str)
        files = contents(~isdir);

        for i = 1:numel(files)
            filename = files(i).name;
            if any(regexp(filename, exclude_dir_with_file_str, 'start'))
              return
            end
        end
    end
    
    % Add d to the path
    p = [d, pathsep];
    
    % Recursively descend through directories
    dirs = contents(isdir); % select only directory entries from the current listing

    for i = 1:numel(dirs)
        dirname = dirs(i).name;
        if ~any(regexp(dirname, exclude_dir_str, 'start'))
          
            % recursive calling of this function.
            p = [p, genpath_exclude_core(fullfile(d, dirname), exclude_dir_str, exclude_dir_with_file_str)];
        end
    end
end
