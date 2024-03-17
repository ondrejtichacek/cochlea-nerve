function [ legendtext ] = Johnson_2011_plot( fig, s, xfun )
%JOHNSON_2011_PLOT 

if nargin < 3
    xfun = @(x) x;
end

legendtext = {};

for k = 1:numel(fig.(s).datasets)

    opt = fig.(s).datasets(k).plotopt;

    x = fig.(s).datasets(k).data(:,1);
    y = fig.(s).datasets(k).data(:,2);

    ncol = size(fig.(s).datasets(k).data,2);

    if ncol == 2
        plot(xfun(x), y, opt{:})
    elseif ncol == 3
        err = fig.(s).datasets(k).data(:,3);

        errorbar(xfun(x), y, err, opt{:})
    elseif ncol == 4
        neg = fig.(s).datasets(k).data(:,3);
        pos = fig.(s).datasets(k).data(:,4);

        errorbar(xfun(x), y, neg, pos, opt{:})
    elseif ncol == 6
        yneg = fig.(s).datasets(k).data(:,3);
        ypos = fig.(s).datasets(k).data(:,4);
        xneg = fig.(s).datasets(k).data(:,5);
        xpos = fig.(s).datasets(k).data(:,6);

        errorbar(xfun(x), y, yneg, ypos, xfun(xneg), xfun(xpos), opt{:})
    else
        error('Number of columns in the dataset = %d is not supported', ncol);
    end

    legendtext{end+1} = fig.(s).datasets(k).name;

end

if numel(legendtext) > 1
    legend(legendtext, 'Location', 'best')
end

xlabel(sprintf('%s [%s]', ...
    fig.(s).xlabel, ...
    fig.(s).xunit))
ylabel(sprintf('%s [%s]', ...
    fig.(s).ylabel, ...
    fig.(s).yunit))

if ~isempty(fig.(s).axis)
    set(gca, fig.(s).axis{:});
end

end

