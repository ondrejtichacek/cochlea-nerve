function [ status, message, messageid ] = srmdir( varargin )
%SREMOVE Safe folder remove

if nargin == 0
    error('Not enough input arguments')
elseif nargin > 2
    error('Too many input arguments')
else
    file = varargin{1};
    if ~isdir(file)
        error('Path %s is not a directory', file)    
    end
    
    file = [file, filesep];
    file = fullfile(file);
    file = replace_tilde(file);
    check_name(file)
    
end

is_in_folders = false;

opt = globalOpt();

folders = [{opt.resDir, opt.tmpDir, opt.scratchdir}, opt.allowedFolders];
for i = 1:length(folders)
    
    folder = fullfile([folders{i}, filesep]);
    folder = replace_tilde(folder);
    
    n = length(folder);
    m = length(file);
    
    if n <= m
        if strcmp(file(1:n), folder)
            is_in_folders = true;
        end
    end
end

if is_in_folders == false
    error('Path %s is not in any of the allowed folders:\n%s', file, ...
        sprintf('%s\n', folders{:}));
else
    if exist(file, 'dir')
        %fprintf('rmdir(%s)\n', file);
        if nargout == 3
            [status, message, messageid] = rmdir(varargin{:});
        else
            rmdir(varargin{:});
        end
    else
        error('Can not remove %s', file)
    end
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

