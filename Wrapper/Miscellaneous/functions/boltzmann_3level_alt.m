function [ y ] = boltzmann_3level_alt( x, M, H1, H2, S1, S2 )

check_size(x, M);
check_size(x, H1);
check_size(x, H2);
check_size(x, S1);
check_size(x, S2);

y = M ./ (1 + 0.5*exp(-(x-H1)./S1) + 0.5*exp(-(x-H2)./S2));

    function check_size(x, V)
        if numel(V) > 1
            assert(all(size(x) == size(V)));
        end
    end

end