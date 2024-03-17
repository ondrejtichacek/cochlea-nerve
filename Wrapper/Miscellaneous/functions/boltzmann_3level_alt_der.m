function [ y ] = boltzmann_3level_alt_der( x, M, H1, H2, S1, S2 )

check_size(x, M);
check_size(x, H1);
check_size(x, H2);
check_size(x, S1);
check_size(x, S2);

y = M * ((exp((H1 - x)/S1))/S1/2 + (exp((H2 - x)/S2))/S2/2) ...
    .* (boltzmann_3level_alt(x, 1, H1, H2, S1, S2).^2);

    function check_size(x, V)
        if numel(V) > 1
            assert(all(size(x) == size(V)));
        end
    end

end