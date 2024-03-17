function [ y ] = boltzmann3_der( x, M, H1, H2, S1, S2 )
%BOLTZMANN_DER

check_size(x, M);
check_size(x, H1);
check_size(x, H2);
check_size(x, S1);
check_size(x, S2);

y = M .* -(-(exp((H1 - x)./S1) .* (exp((H2 - x)./S2) + 1))./S1 - exp((H1 - x)./S1 + (H2 - x)./S2)./S2) ./ (exp((H1 - x)./S1).*(exp((H2 - x)./S2) + 1) + 1).^2;

    function check_size(x, V)
        if numel(V) > 1
            assert(all(size(x) == size(V)));
        end
    end

end

