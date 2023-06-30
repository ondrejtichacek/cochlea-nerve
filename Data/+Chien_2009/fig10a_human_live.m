function [xnew, y, yerr, ymag, yerrmag] = fig10a_human_live(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(300),log10(6000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

f_mean = fullfile(folder, 'f10a-human-live.csv');
f_top_err = fullfile(folder, 'f10a-human-live-top_err.csv');
f_bot_err = fullfile(folder, 'f10a-human-live-bot_err.csv');


data_x_lim = [300, 6000];
if any(xnew < data_x_lim(1)) || any(xnew > data_x_lim(2))
    disp(data_x_lim);
    disp([min(xnew), max(xnew)]);
    warning('Frequencies out of range!')
end

tmp = readmatrix(f_mean);
x = tmp(:,1);
y = tmp(:,2);

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    plot(x, y, 'kx-');
    ylabel('magnitude [\mu m /s /Pa]')
    xlabel('frequency (Hz)')
end
    
Fmean = griddedInterpolant(x, y, 'linear', 'linear');

tmp = readmatrix(f_top_err);
xet = tmp(:,1);
yet = tmp(:,2);

if args.plotflag
    plot(xet, yet, 'kx-');
end

Ferrtop = griddedInterpolant(xet, yet, 'linear', 'linear');

tmp = readmatrix(f_bot_err);
xeb = tmp(:,1);
yeb = tmp(:,2);

if args.plotflag
    plot(xeb, yeb, 'kx-');
end

Ferrbot = griddedInterpolant(xeb, yeb, 'linear', 'linear');

y = Fmean(xnew);

yet = -y + Ferrtop(xnew);
yeb = y - Ferrbot(xnew);

yerr = [yet; yeb];

if args.plotflag
    plot(xnew, Fmean(xnew), 'bo', ...
         xnew, Ferrtop(xnew), 'ro', ...
         xnew, Ferrbot(xnew), 'ro');
end


mag2db = @(z) 20*log10(z);

yerrmag = [ ...
            mag2db(y + yerr(1,:)) - mag2db(y); ...
            mag2db(y) - mag2db(y - yerr(2,:))];
        
ymag = mag2db(y);
ymag = ymag - max(ymag);

end

