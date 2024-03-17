function p = genpath_exclude(d, exclude_dir_pattern, exclude_dir_with_file_pattern)
arguments
    d (1,:) char
    exclude_dir_pattern cell {iscellstr} = {}
    exclude_dir_with_file_pattern cell {iscellstr} = {}
end

    if ~isempty(exclude_dir_pattern)
        default_exclude_pattern = { ...
            '\.', ...
            '\.\.', ...
            'private', ...
            '\@.*', ...
            '\+.*', ...
            'resources' };

        exclude_dir_pattern = [exclude_dir_pattern, default_exclude_pattern];

        for i = 1:numel(exclude_dir_pattern)
            exclude_dir_pattern{i} = ['^', exclude_dir_pattern{i}, '$'];
        end
        
        exclude_dir_str = join(exclude_dir_pattern, '|');
        exclude_dir_str = exclude_dir_str{1};
    else
        exclude_dir_str = '';
    end
    
    if ~isempty(exclude_dir_with_file_pattern)
        for i = 1:numel(exclude_dir_with_file_pattern)
            exclude_dir_with_file_pattern{i} = ['^', exclude_dir_with_file_pattern{i}, '$'];
        end
        
        exclude_dir_with_file_str = join(exclude_dir_with_file_pattern, '|');
        exclude_dir_with_file_str = exclude_dir_with_file_str{1};
    else
        exclude_dir_with_file_str = '';
    end
    
    p = genpath_exclude_core(d, exclude_dir_str, exclude_dir_with_file_str);

end
