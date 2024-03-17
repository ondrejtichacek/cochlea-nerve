function x = gaussgrid(N, plotflag)
%GAUSSGRID Generates a set of points covering the interval 0-1 with Gaussian distributed density
%  N = Number of points in the grid
%  X = BM site array of gaussian distributed spacing
%  PLOT_FLAG = plot only if 1
%
arguments
    N = 500;
    plotflag = 0;
end

stred = 0.35;
sig1 = 2*(1.2)^2;
sig2 = 2*(0.46)^2;

xo = (1:N)/N;

ind = ceil(stred/xo(1));
dx(1:ind) = exp(((xo(1:ind)-xo(ind)).^2)./sig1);
dx(ind+1:N) = exp(((xo(ind+1:N)-xo(ind)).^2)./sig2);
dx = dx/sum(dx);

x(1) = dx(1);

for i = 2:N
    x(i) = x(i-1) + dx(i);
end

if x(N) > 1
    x(N) = 1;
end

dens = ones(size(dx))./dx;
Nmx = max(dens);
Nmn = min(dens);

if plotflag == 1
    [~, im] = min(abs(x-0.5));
    Nmdl = dens(im);

    figure;
    
    plot(x, zeros(size(x)) + Nmx/2,'.w',x,dens,'-r'),
%     axis([0,1,0,Nmx]);
    text(0.2, Nmn + 0.45*(Nmx-Nmn), sprintf('Max. density = %.4g', Nmx));
    text(0.2, Nmn + 0.35*(Nmx-Nmn), sprintf('Mid. density = %.4g', Nmdl));
    text(0.2, Nmn + 0.25*(Nmx-Nmn), sprintf('Min. density = %.4g', Nmn));
    ylabel('Point density')
    xlabel('Fractional distance from stapes')
end

end