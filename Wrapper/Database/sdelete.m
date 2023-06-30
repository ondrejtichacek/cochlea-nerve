function sdelete( varargin )
%SREMOVE Safe file remove

if nargin == 0
    error('Not enough input arguments')
end

files = varargin;
is_in_folders = logical(zeros(size(files)));

for j = 1:nargin
    file = fullfile(files{j});
    
    if ~exist(file, 'file')
        if exist(file, 'dir')
            error('Path %s is not a file', file);
        else
            error('Path %s does not exist', file);
        end
    end
    file = replace_tilde(file);
    check_name(file)
    files{j} = file;
end

opt = globalOpt();

folders = [{opt.resDir, opt.tmpDir, opt.scratchdir}, opt.allowedFolders];
for i = 1:length(folders)
    
    folder = fullfile([folders{i}, filesep]);
    folder = replace_tilde(folder);
    
    n = length(folder);
    for j = 1:nargin
        file = files{j};
        m = length(file);

        if n <= m
            if strcmp(file(1:n), folder)
                is_in_folders(j) = true;
            end
        end
    end
end

if any(is_in_folders == false)
    for j = 1:nargin
        if is_in_folders(j) == false
            error('Path %s is not in any of the allowed folders:\n%s', files{j}, ...
                sprintf('%s\n', folders{:}));
        end
    end
end

if all(is_in_folders)
    % fprintf('delete(%s)\n', varargin{:});
    delete(varargin{:})
end

function check_name(folder)
    if any(strfind(folder,'..'))
        error('Path %s contains `..`, which is not supported', folder);
    end
end
    
function folder = replace_tilde(folder)
    if find(folder == '~', 1)
        if isunix()
            [st, full_home_path] = system('echo $HOME');
            if st ~= 0
                error('Error when executing `echo $HOME`');
            end
            full_home_path = full_home_path(1:end-1);
            folder = strrep(folder, '~', full_home_path);
        else
            error('Folder path contains `~`');
        end
    end
end

end

