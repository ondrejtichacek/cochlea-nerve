function [x, y, xerr, yerr] = combineDatasets(datasets, varargin)

    specific_datasets = varargin;

    if isempty(specific_datasets)
        data = [datasets.data];
    else
        I = zeros(1, numel(specific_datasets));
        for k = 1:numel(specific_datasets)
            I(k) = find(strcmp(specific_datasets{k}, {datasets.name}));
        end
        data = vertcat(datasets(I).data);
    end

    [x, ind] = sort(data(:,1));
    y = data(ind,2);

    if nargout > 2

        switch size(data,2)
            case 3 % x, y, yerr
                y_minus = y - data(ind,3);
                y_plus = y + data(ind,3);
            case 4 % x, y, y_err_min, y_err_max
                y_minus = y - data(ind,3);
                y_plus = y + data(ind,4);
        end

        xerr = [x, x];
        yerr = [y_minus, y_plus];
    end


end