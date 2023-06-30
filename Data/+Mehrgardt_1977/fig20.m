function [xnew, ya, yb] = fig20(xnew, args)
arguments
    xnew (1,:) double = logspace(log10(200),log10(15000),200);
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

f_a = fullfile(folder, 'f20a.csv');
f_b = fullfile(folder, 'f20b.csv');


data_x_lim = [200, 15000];
if any(xnew < data_x_lim(1)) || any(xnew > data_x_lim(2))
    disp(data_x_lim);
    disp([min(xnew), max(xnew)]);
    warning('Frequencies out of range!')
end

tmp = readmatrix(f_a);
xa = tmp(:,1);
ya = tmp(:,2);

tmp = readmatrix(f_b);
xb = tmp(:,1);
yb = tmp(:,2);


if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    plot(xa, ya);
    plot(xb, yb);
end

Fa = griddedInterpolant(xa, ya, 'linear', 'none');
Fb = griddedInterpolant(xb, yb, 'linear', 'none');


% offset = 1.5;
offset = 0;

ya = Fa(xnew) + offset;
yb = Fb(xnew) + offset;

end

