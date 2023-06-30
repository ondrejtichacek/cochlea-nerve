function [ tspan_required ] = findTspan( PART, resultFound, resultPath, topt, args )
%FINDTSPAN
arguments
    PART (1,:) char
    resultFound (1,1) logical
    resultPath (1,:) char
    topt (1,1)
    args.midopt = []
    args.mechopt = []
    args.mnaopt = []
end

time_tolerance = Time(1, 'ns');

if ~resultFound
    tspan_required = topt;
    
else
    found = load(fullfile(resultPath, 'settings.mat'), 'topt');

    tf_required = topt.total.tf;
    tf_found  = found.topt.total.tf;

    if tf_required <= tf_found
        % no need to compute anything
        tspan_required = [];
    else
        % need to compute the result from `tf_found` to `tf_required`
        tspan_required = tspanOpt('t0', tf_found, 'tf', tf_required);
    end

    switch PART
        case 'middleear'
            assert(isa(args.midopt, 'midOpt'))
            
        case 'mechanical'
            assert(isa(args.mechopt, 'mechOpt'))
            
        case 'electrical'
            
            assert(isa(args.mnaopt, 'mnaOpt'))
            
            time_file = fullfile(resultPath, 'Time.mat');
            if exist(time_file, 'file')
                mf = matfile(time_file, 'Writable', false);
                
%                 if any([mf.T(1,1) ~= found.topt.total.t0, ...
%                         mf.T(end,1) ~= found.topt.total.tf])
%                     fprintf('Deleting %s\n', resultPath);
%                     deleteResultFiles(resultPath);
%                 end
                
                assert(isequaleps(time_tolerance, Time(mf.T(1,1),mf.unit), found.topt.total.t0), 'Time points do not correspond.');
                assert(isequaleps(time_tolerance, Time(mf.T(end,1),mf.unit), found.topt.total.tf), 'Time points do not correspond.');
                T_size = max(matfileWhosByName(time_file, 'T', 'size'));
            end

            val_files = {fullfile(resultPath, 'Volt.mat'), fullfile(resultPath, 'Curr.mat')};
            for i = 1:numel(val_files)
                if exist(val_files{i}, 'file')
                    var_sizes = matfileWhosByName(val_files{i}, [], 'size');
                    vars = fields(var_sizes);
                    for j = 1:numel(vars)
                        assert(all(sort(var_sizes.(vars{j})) == sort([T_size, args.mnaopt.Numstacks])), ...
                            sprintf('Variable %s in the file %s has a wrong shape', vars{j}, val_files{i}));
                    end
                end
            end
        case 'synapse'
        case 'nerve'
        otherwise
            error('PART %s is not recognised option.',PART);
    end
end

end

