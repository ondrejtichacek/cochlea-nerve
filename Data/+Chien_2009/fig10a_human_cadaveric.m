function [xnew, y, ymag] = fig10a_human_cadaveric(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(50),log10(10000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

f_mean = fullfile(folder, 'f10a-human-cadaveric.csv');

data_x_lim = [50, 10000];
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

y = Fmean(xnew);

mag2db = @(z) 20*log10(z);

ymag = mag2db(y);
ymag = ymag - max(ymag);

end

