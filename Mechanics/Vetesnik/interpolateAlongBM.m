function [ fun ] = interpolateAlongBM(lambda, opts)
%INTERPOLATEALONGBM
arguments
    lambda
    opts.x (1,:) double = []
    opts.cf (1,:) double = []
    opts.plotflag (1,1) logical = false
    opts.hfig = []
end

if opts.plotflag == true && isempty(opts.hfig)
    opts.hfig = figure;
end

if isempty(opts.cf) && isempty(opts.x)
    error('Setting either cf or x is required');
end
if ~isempty(opts.cf) && ~isempty(opts.x)
    error('It is possible to set either cf or x, not both');
end

if isempty(opts.x)
    [cf, i] = sort(opts.cf, 'descend');
    lambda = lambda(i);

    cf_model = 'loglinear';
    cf_parameters = [4.3647; -2.8577];
    
    z = characteristic_position(cf, cf_model, cf_parameters);
else
    [z, i] = sort(opts.x, 'ascend');
    lambda = lambda(i);
end

% L = griddedInterpolant(z, lambda, 'pchip', 'linear');
L = griddedInterpolant(z, lambda, 'pchip');

fun = @(x) L(x);

if opts.plotflag == true
    actfigure(opts.hfig);
    x = linspace(0,1,300);
    hold on
    plot(z, lambda, 'bo')
    plot(x,fun(x))
    xlabel('norm. distance from stapes')
end

end
