function dispstruct(S)

T = S;

f = fieldnames(T);

for i = 1:numel(f)
    if isa(T.(f{i}), 'Unit')
        if numel(T.(f{i})) == 0
            T.(f{i}) = sprintf( ...
                '[]   [0x0 %s]', ...
                class(T.(f{i})));        
        elseif numel(T.(f{i})) == 1
            T.(f{i}) = sprintf( ...
                '%s   [1x1 %s]', ...
                T.(f{i}), class(T.(f{i})));        
        end
    end
end

disp(T);

end

