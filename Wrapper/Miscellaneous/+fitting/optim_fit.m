function [fun, xbest, fval] = optim_fit(FUN, tdata, ydata, x0, opts, tfun, do_fit)

if nargin < 4 || isempty(x0)
    x0 = rand(4,1);
else
    x0(isnan(x0)) = rand(sum(isnan(x0)),1);
end

if nargin < 5 || isempty(opts)
    opts = optimset('MaxFunEvals', 1e3*numel(x0), 'MaxIter', 1e3*numel(x0));
end

if nargin < 6 || isempty(tfun)
    tfun = @(t) t;
end

if nargin < 7 || isempty(do_fit)
    do_fit = true;
end

%fitfun = @(x,t) x(1) + FUN(tfun(t), x(2), x(3), x(4));

if do_fit == true
    objfun = @(x) sum((ydata - fitfun(x, tdata, FUN, tfun)).^2);

    [xbest, fval] = fminsearch(objfun, x0, opts);
    
else    
    xbest = x0;
    fval = NaN;
end

fun = @(t) fitfun(xbest, t, FUN, tfun);

end

% this needs to be outside of the workspace of optim_fit otherwise we have
% a problem in hashing ...
function F = fitfun(x,t, FUN, tfun)
xx = num2cell(x(2:end));
F = x(1) + FUN(tfun(t), xx{:});
end