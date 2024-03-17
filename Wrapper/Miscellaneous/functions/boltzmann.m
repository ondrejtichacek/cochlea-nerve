function [ y ] = boltzmann( x, M, H, S )
%BOLTZMANN 

check_size(x, M);
check_size(x, H);
check_size(x, S);

y = M./(1+exp(-(x-H)./S));

    function check_size(x, V)
        if numel(V) > 1
            assert(all(size(x) == size(V)));
        end
    end

end

