function [xnew, y, yerr] = fig3a(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(100),log10(20000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

%% Magnitude

f_mean = fullfile(folder, 'f3a_raw_mean.csv');
f_top_err = fullfile(folder, 'f3a_raw_top_err.csv');

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


tmp = readmatrix(f_top_err);
x = tmp(:,1);
y = tmp(:,2);

if args.plotflag
    plot(x, y);
end

Ferr = griddedInterpolant(x, y);

y = Fmean(xnew);
yerr = Ferr(xnew) - y;


end

