function erbticks(I)
arguments
    I = 1:9
end

t = 2.^(I);

b = 125;

set(gca, ...
    'Xtick', b*t, ...
    'XTickLabels', cellfun(@(x) sprintf('%d', round(x)), num2cell(b*t), 'UniformOutput', false));

end