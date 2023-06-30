function [ s, x ] = slicesToExcite( slices, Numstacks, Volt, results, xgrid )
%SLICESTOEXCITE parse options to determine which slices to excite

x = xgrid;

if ~isempty(results)
    jj = find(strcmp({results.analysis.fullid}, 'IHCabsolute'), 1);
    if isempty(jj)
        error('Expected to find a file with variable fullid = %s', 'IHCabsolute')
    end
    % stats = load(results.analysis(jj).analysis_file);
    stats = matfile(results.analysis(jj).analysis_file, 'Writable', false);
    % x = results.result_files.fBM.x;
else
    stats = [];
    % x = [];
end

if length(slices) == 1
    if ischar(slices{1})
        switch slices{1}
            case 'all'
                s = 1:Numstacks;
            case 'max'
                [~,ii_max] = max(max(abs(Volt - Volt(1,:)),[],1));
                s = ii_max;
            case 'max_std'
                [~,ii_max] = max(std(Volt,[],1));
                s = ii_max;
            case 'max_std_new'
                tmp = stats.oscillation_max;
                s = tmp.spectrum_restricted.ind;
            otherwise
                error('Unsupported option');
        end  
    elseif isnumeric(slices{1})
        s = slices{1};
    else
        error('Unsupported option type')
    end
elseif length(slices) == 2
    if ischar(slices{1}) && isnumeric(slices{2})
        
        [~,ii_max] = max(max(Volt,[],1));
        
        switch slices{1}
            case 'max+'
                s = ii_max:ii_max+slices{2};
            case 'max-'
                s = ii_max-slices{2}:ii_max;
            case 'max+-'
                s = ii_max-slices{2}:ii_max+slices{2};
            otherwise
                error('Unsupported option');
        end  
    else
        error('Unsupported option type')
    end
end

if s(1) < 1
    s = 1:s(end);
    warning('Slice sellection out of range, setting slicesToExcite = 1 : %d', s(end));    
end
if s(end) > Numstacks
    s = s(1):Numstacks;
    warning('Slice sellection out of range, setting slicesToExcite = %d : %d', s(1), Numstacks);    
end

if ~isempty(x)
    x = x(s);
end

end

