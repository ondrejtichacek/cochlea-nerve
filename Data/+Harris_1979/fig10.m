function [x, y, z] = fig10(args)
%FIG10 
arguments
    args.plotflag (1,1) logical = false
    args.plot_fit (1,1) logical = true
end

p10_1 = [0.96980, 80.770
1.9356, 82.103
4.8609, 82.194
9.7107, 86.228
19.647, 89.984
48.686, 97.555
97.219, 100.34
NaN, NaN
NaN, NaN
];

p10_2 = [0.94267, 80.816
1.8895, 82.123
4.7181, 82.157
9.3951, 86.169
18.955, 89.967
47.023, 97.617
93.636, 100.42
NaN, NaN
NaN, NaN
];

p20_1 = [0.99325, 48.707
1.9701, 50.524
4.9071, 54.563
9.9002, 59.358
19.702, 70.594
49.434, 86.336
97.433, 97.502
196.48, 100.91
NaN, NaN
];

p20_2 = [0.97091, 48.717
1.9334, 50.523
4.7961, 54.756
9.6765, 59.551
19.143, 70.751
47.800, 86.372
94.252, 97.643
189.54, 101.09
NaN, NaN
];

p30_1 = [0.98681, 30.216
1.9693, 31.203
4.8913, 36.280
9.8355, 40.729
19.564, 50.718
48.684, 70.131
96.922, 82.544
194.72, 93.572
290.55, 95.174
];

p30_2 = [0.97410, 30.141
1.9525, 31.235
4.8277, 36.322
9.6765, 40.690
19.269, 50.751
47.332, 70.358
94.252, 82.625
187.68, 93.612
281.88, 95.122
];

p40_1 = [0.98634, 19.759
1.9659, 26.355
4.8918, 27.485
9.8006, 30.618
19.447, 42.685
48.905, 64.868
97.675, 77.281
195.75, 90.386
292.20, 93.028
];

p40_2 = [0.98051, 19.821
1.9525, 26.466
4.8435, 27.568
9.7083, 30.654
19.143, 42.637
47.644, 64.877
94.872, 77.216
188.92, 90.552
281.88, 93.129
];

p60_1 = [0.98793, 15.258
1.9715, 16.244
4.8767, 18.690
9.7868, 26.602
19.637, 33.959
48.722, 54.203
98.166, 64.262
195.00, 79.445
292.06, 82.640
];

p60_2 = [0.98051, 15.337
1.9589, 16.359
4.8754, 18.742
9.7083, 26.668
19.396, 34.025
47.800, 54.201
95.810, 64.405
189.54, 79.449
283.73, 82.667
];


% data from table 1
r_d = [22, 53, 71, 81, 86];
gamma = [32, 37, 55, 66, 74];


override_x = [1,2,5,10,20,50,100,200,300];

z = [10, 20, 30, 40, 60];

p10 = mean(cat(3, p10_1, p10_2),3);
p20 = mean(cat(3, p20_1, p20_2),3);
p30 = mean(cat(3, p30_1, p30_2),3);
p40 = mean(cat(3, p40_1, p40_2),3);
p60 = mean(cat(3, p60_1, p60_2),3);


y = [p10(:,2), p20(:,2), p30(:,2), p40(:,2), p60(:,2)];

x = override_x;

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    
    cmap = lcmp(numel(z));
    xx = logspace(log10(1), log10(300), 100);

    for i = 1:numel(z)
        l = z(i);
        
        if args.plot_fit
            ls = 'o';
        else
            ls = 'o:';
        end

        plot(x, y(:,i), ls, ...
            'Color', cmap(i,:), ...
            'DisplayName', sprintf('+%d dB', l))
    
        if args.plot_fit
            plot(xx, 100-r_d(i)*exp(-xx/gamma(i)), '-', ...
                'Color', cmap(i,:), ...
                'DisplayName', '')
        end
    end
    xlabel('\Delta T (ms)')
    ylabel('Median probe response magnitude (% control)')
    legend('Location', 'best')

    xlim([0.9,330])
    ylim([10, 105])
end

end

