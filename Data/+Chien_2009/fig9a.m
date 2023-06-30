function [xnew, y, yerr, ymag, yerrmag] = fig9a(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(100),log10(6000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

f_mean = fullfile(folder, 'f9a_cadaver_raw_mean.csv');
f_top_err = fullfile(folder, 'f9a_cadaver_raw_top_err.csv');
f_bot_err = fullfile(folder, 'f9a_cadaver_raw_bottom_err.csv');


data_x_lim = [100, 6000];
if any(xnew < data_x_lim(1)) || any(xnew > data_x_lim(2))
    disp(data_x_lim);
    disp([min(xnew), max(xnew)]);
    warning('Frequencies out of range, using nearest (constant) extrapolation')
end

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
    
Fmean = griddedInterpolant(x, y, 'linear', 'nearest');

tmp = readmatrix(f_top_err);
xet = tmp(:,1);
yet = tmp(:,2);

if args.plotflag
    plot(xet, yet, '.');
end

Ferrtop = griddedInterpolant(xet, yet, 'linear', 'nearest');

tmp = readmatrix(f_bot_err);
xeb = tmp(:,1);
yeb = tmp(:,2);

if args.plotflag
    plot(xeb, yeb, '.');
end

Ferrbot = griddedInterpolant(xeb, yeb, 'linear', 'nearest');

y = Fmean(xnew);

yet = -y + Ferrtop(xnew);
yeb = y - Ferrbot(xnew);

yerr = [yet; yeb];

if args.plotflag
    plot(xnew, Ferrtop(xnew), xnew, Ferrbot(xnew), '--');
end


mag2db = @(z) 20*log10(z);

yerrmag = [ ...
            mag2db(y + yerr(1,:)) - mag2db(y); ...
            mag2db(y) - mag2db(y - yerr(2,:))];
        
ymag = mag2db(y);
ymag = ymag - max(ymag);

end

