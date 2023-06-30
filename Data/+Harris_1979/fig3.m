function [x, y] = fig3(args)
%FIG3 
arguments
    args.plotflag (1,1) logical = false
end

folder = fileparts(mfilename('fullpath'));

f = fullfile(folder, 'fig3.csv');


tmp = readmatrix(f);
x = tmp(:,1);
y = tmp(:,2);

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    plot(x, y, 'bo-');
    xlabel('\Delta T (ms)')
    ylabel('Probe response magnitude (% control)')
end

end

