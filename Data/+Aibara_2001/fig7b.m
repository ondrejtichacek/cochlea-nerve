function [xnew, y, yerr] = fig7b(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(50),log10(10000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

f_mean = fullfile(folder, 'f7b_raw_mean.csv');
f_bottom_err = fullfile(folder, 'f7b_raw_bottom_err.csv');


tmp = readmatrix(f_mean);
x = tmp(:,1);
y = tmp(:,2);

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    plot(x, y);
end
    
Fmean = griddedInterpolant(x, y);


tmp = readmatrix(f_bottom_err);
x = tmp(:,1);
y = tmp(:,2);

if args.plotflag
    plot(x, y);
end

Ferr = griddedInterpolant(x, y);

y = Fmean(xnew);

yerr = y - Ferr(xnew);



end

