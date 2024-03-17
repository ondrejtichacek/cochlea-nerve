function [ undamp ] = undamping( bigamma, damp0, gain, lambda, x )
%UNDAMPING
arguments
    bigamma (:,1) double
    damp0 (:,1) double
    gain (1,1) double
    lambda
    x  (:,1) double
end


switch class(lambda)
    case 'char'
        mf = matfile(lambda);
        F = griddedInterpolant(mf.x, mf.lam);
        lam = F(x);
    case 'double'
        lam = lambda * ones(size(x));
    case 'function_handle'
        lam = lambda(x);
    otherwise
        error('Parameter lambda must be either char, double or function_handle.')
end

undamp = gain*damp0.*bigamma.*lam;

end

