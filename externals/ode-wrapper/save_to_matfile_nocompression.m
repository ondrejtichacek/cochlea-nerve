function save_to_matfile_nocompression(mfile, varname, variable, args)
%SAVE_TO_MATFILE_NOCOMPRESSION
arguments
    mfile
    varname (1,:) char
    variable
    args.concatenate (1,1) logical = true
    args.concat_dim (1,1) double = 1
end

saved = false;

file_exists = exist(mfile.Properties.Source, 'file');
            
if ~file_exists || isempty(who(mfile, varname))

    if file_exists
        saveflags = {'-append'}; % v7.3 not included to prevent warning
    else
        saveflags = {'-v7.3'};
    end

    VAR.(varname) = variable;

    save(mfile.Properties.Source, ...
        '-struct', 'VAR', ...
        '-nocompression', saveflags{:});

    % mfile = matfile(mfile.Properties.Source, ...
    %     'Writable', true);

    saved = true;
end

if ~saved
    
    assert(args.concatenate == true)
                
    s = matfileWhosByName(mfile, varname, 'size');

    assert(numel(s) <= 2); % only for 1-2 dim variables
    
    n0 = s(args.concat_dim);
    
    assert(n0 > 0)

    int = 1:size(variable, args.concat_dim);

    if args.concat_dim == 1
        mfile.(varname)(n0 + int, :) = variable;
        
    elseif args.concat_dim == 2
        mfile.(varname)(:, n0 + int) = variable;
        
    else
        error('not implemented')
    end
end

end

