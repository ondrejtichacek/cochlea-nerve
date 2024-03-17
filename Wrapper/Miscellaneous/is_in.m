function [TF] = is_in(val, cellset)

    TF = false;

    for i = 1:numel(cellset)
        if isstring(val) && isstring(cellset{i})
            if string(val) == string(cellset{i})
                TF = true;
                return
            end
        else
            if all(size(val) == size(cellset{i}))
                if val == cellset{i}
                    TF = true;
                    return
                end
            end
        end
    end
end