function [x, y, z] = fig15a(args)
%FIG10 
arguments
    args.plotflag (1,1) logical = false
end

d2 = [
    1.0138, 87.199
2.0290, 88.268
5.0171, 90.184
10.097, 94.424
20.220, 97.103
50.450, 100.15
NaN, NaN
NaN, NaN
];

d5 = [
    1.0188, 76.978
2.0312, 79.183
5.0239, 81.998
10.114, 87.232
20.197, 93.932
50.427, 98.924
NaN, NaN
NaN, NaN
];

d10 = [
    1.0077, 70.401
2.0087, 72.133
5.0394, 78.544
10.183, 82.074
20.166, 89.579
50.630, 98.119
101.37, 99.993
NaN, NaN
];

d50 = [
    1.0025, 56.016
2.0094, 61.108
5.0367, 65.058
10.156, 74.881
20.206, 83.144
50.573, 94.949
101.37, 99.993
NaN, NaN
];

d100 = [
    1.0125, 47.877
2.0012, 49.799
5.0167, 53.986
10.142, 58.935
20.284, 69.894
50.397, 85.249
100.83, 97.202
202.80, 99.927
];

d200 = [
    1.0129, 48.871
2.0106, 50.840
5.0402, 54.979
10.146, 60.071
20.126, 72.071
50.433, 87.236
101.18, 94.930
202.80, 99.927
];

% data from table 2
r_d = [14, 25, 31, 44, 53, 53];
gamma = [12, 14, 18, 23, 37, 37];


override_x = [1,2,5,10,20,50,100,200];

z = [2, 5, 10, 50, 100, 200];


y = [d2(:,2), d5(:,2), d10(:,2), d50(:,2), d100(:,2), d200(:,2)];

x = override_x;

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    
    cmap = lcmp(numel(z));
    xx = logspace(log10(1), log10(200), 100);

    for i = 1:numel(z)
        l = z(i);
        
        plot(x, y(:,i), 'o', ...
            'Color', cmap(i,:), ...
            'DisplayName', sprintf('+%d dB', l))
    
        plot(xx, 100-r_d(i)*exp(-xx/gamma(i)), '-', ...
            'Color', cmap(i,:), ...
            'DisplayName', '')
    end
    xlabel('\Delta T (ms)')
    ylabel('Median probe response magnitude (% control)')
    legend('Location', 'best')

    xlim([0.9,330])
    ylim([10, 105])
end

end

