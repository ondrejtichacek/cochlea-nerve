function [xnew, y, yerr] = fig9b(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(100),log10(6000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

%% Phase lag

f_mean = fullfile(folder, 'f9b_cadaver_raw_mean.csv');
f_top_err = fullfile(folder, 'f9b_cadaver_raw_top_err.csv');
f_bot_err = fullfile(folder, 'f9b_cadaver_raw_bottom_err.csv');

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
xet = tmp(:,1);
yet = tmp(:,2);

if args.plotflag
    plot(xet, yet, '.');
end

tmp = readmatrix(f_bot_err);
xeb = tmp(:,1);
yeb = tmp(:,2);

if args.plotflag
    plot(xeb, yeb, '.');
end

yet = yet - Fmean(xet);
yeb = yeb - Fmean(xeb);

y = [-yeb; yet];
x = [xeb; xet];

[x, i] = sort(x);
y = y(i);

Ferr = griddedInterpolant(x, y);

y = Fmean(xnew);
yerr = Ferr(xnew);

if args.plotflag
    plot(xnew, y - yerr, xnew, y + yerr, '--');
end

end