function [c, cmap] = linecolors(i)
%LINECOLORS 

colors =  [ % matlab default colors
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

nc = size(colors, 1);

c = colors(1 + mod(i-1,nc),:);

cmap = repmat(colors,ceil(i/nc), 1);
cmap = cmap(1:i, :);

end

