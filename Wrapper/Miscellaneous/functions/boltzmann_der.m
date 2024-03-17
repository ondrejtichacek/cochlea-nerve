function [ y ] = boltzmann_der( x, M, H, S )
%BOLTZMANN_DER

check_size(x, M);
check_size(x, H);
check_size(x, S);

y = M .* exp((H + x)./S) ./ (S.*(exp(H./S) + exp(x./S)).^2);

% alternatively

% y = M./(2*S + 2*S.*cosh((H - x)./S));

    function check_size(x, V)
        if numel(V) > 1
            assert(all(size(x) == size(V)));
        end
    end

end

