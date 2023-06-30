function [xnew, y, yerr, ymag, yerrmag] = fig7a(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(50),log10(10000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

f_mean = fullfile(folder, 'f7_raw_mean.csv');
f_top_err = fullfile(folder, 'f7_raw_top_err.csv');


tmp = readmatrix(f_mean);
x = tmp(:,1);
y = tmp(:,2);

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
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

ytop = Ferr(xnew);
yerr = ytop - y;

mag2db = @(z) 20*log10(z);

yerrmag = [ ...
    mag2db(y + yerr) - mag2db(y);
    mag2db(y) - mag2db(y - yerr)];

ymag = mag2db(y);
ymag = ymag - max(ymag);


end