function [ cleanup ] = handleLockFile( path, args )
%HANDLELOCKFILE
arguments
    path char
    args.action string = "check"
end

f_lock = fullfile(path, 'LOCK');


if args.action == "purge"
    
    if exist(f_lock, 'file')
        warning('Deleting lock file');
        delete(f_lock);
    end
elseif any(args.action == ["check", "create"])

    if exist(f_lock, 'file')
        error('Lock file exists: %s\n ==> Aborting.\n', f_lock)
    end
else
    error('Unrecognised action %s', args.action);
end

if args.action == "create"
    if ~exist(path, 'dir')
        error('Can''t create LOCK file - path %s does not exist\n', path)
    end

    % create lock file
    fid = fopen(f_lock, 'w');
    fclose(fid);

    % delete the file on cleanup
    cleanup = onCleanup(@() sdelete(f_lock));
end


end

