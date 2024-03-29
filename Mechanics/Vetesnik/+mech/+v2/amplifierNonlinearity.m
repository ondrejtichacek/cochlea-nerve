function [ RHS ] = amplifierNonlinearity( ~, y, N, ~, TMa, BMy, ~, nonlinearity )

fun  = BMy .* nonlinearity(y(2*N+1:3*N));

RHS = [ ...
    zeros(N,1); ...
    -fun; ...
    zeros(N,1); ...
    TMa.*fun];

end

