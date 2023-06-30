function [x, y, z] = fig14b(args)
%FIG10 
arguments
    args.plotflag (1,1) logical = false
end

d10 = [
1.0004, 84.894
2.0027, 83.188
4.9923, 82.127
9.9681, 84.154
19.938, 95.740
49.524, 99.687
99.100, 98.618
200.12, 96.639
];

d20 = [
1.0030, 57.216
2.0043, 58.059
5.0216, 61.824
9.9491, 74.867
19.945, 83.176
49.452, 89.763
99.030, 99.620
NaN, NaN
];

d50 = [
0.99890, 51.025
1.9971, 51.230
4.9895, 58.910
9.9117, 68.220
19.988, 80.172
49.712, 82.297
99.463, 93.429
199.86, 98.460
];

d100 = [
1.0024, 46.108
2.0008, 48.590
4.9915, 58.364
9.9276, 65.944
19.905, 74.071
49.655, 83.936
98.323, 97.798
199.71, 99.553
];

override_x = [1,2,5,10,20,50,100,200];

z = [10, 20, 50, 100];

y = [d10(:,2), d20(:,2), d50(:,2), d100(:,2)];

x = override_x;

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    
    cmap = lcmp(numel(z));

    for i = 1:numel(z)
        l = z(i);
        
        plot(x, y(:,i), 'o-', ...
            'Color', cmap(i,:), ...
            'DisplayName', sprintf('masker %d ms', l))
   
    end
    xlabel('\Delta T (ms)')
    ylabel('Median probe response magnitude (% control)')
    legend('Location', 'best')

    %xlim([0.9,330])
    %ylim([10, 105])
end

end

