function [s] = update_struct(s, recursive, f, v)
%UPDATE_STRUCT


arguments
    s (1,1) struct
    recursive (1,1) logical
end
arguments (Repeating)
    f (1,:) char
    v
end

for i = 1:numel(f)
    if recursive == true && isa(v{i}, 'struct')
        s.(f{i}) = update_struct(v{i}, recursive, v{i});
    else
        s.(f{i}) = v{i};
    end
end

end

