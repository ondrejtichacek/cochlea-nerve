function [S] = nanstruct(S)
%NANSTRUCT Creates a f
arguments
    S (1,1) struct
end

f = fieldnames(S);
for i = 1:length(f)
    if isnumeric(S.(f{i}))
        S.(f{i}) = NaN(size(S.(f{i})), 'like', S.(f{i}));
    elseif isa(S.(f{i}), 'Unit')
        S.(f{i}) = NaN * S.(f{i});
    elseif isstruct(S.(f{i}))
        S.(f{i}) = nanstruct(S.(f{i}));
    elseif islogical(S.(f{i}))
        S.(f{i}) = NaN * S.(f{i});
    elseif ischar(S.(f{i}))
        S.(f{i}) = S.(f{i});
    else
        error('Not yet implemented for class %s of the structure field %s', class(S.(f{i})), f{i})
    end
end

