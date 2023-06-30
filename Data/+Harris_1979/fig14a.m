function [x, y, z] = fig14a(args)
%FIG10 
arguments
    args.plotflag (1,1) logical = false
end

d2 = [
1.0057, 84.124
1.9828, 83.342
4.9832, 90.064
10.045, 97.046
20.109, 100.10
NaN, NaN
NaN, NaN
NaN, NaN
];

d5 = [
1.0027, 71.155
1.9747, 73.022
5.0018, 81.297
10.074, 90.105
20.117, 99.188
NaN, NaN
NaN, NaN
NaN, NaN
];

d10 = [
1.0059, 63.758
1.9789, 67.908
5.0214, 72.073
10.027, 81.155
20.126, 98.092
50.259, 99.973
NaN, NaN
NaN, NaN
];

d50 = [
1.0095, 55.265
1.9849, 60.876
4.9944, 64.767
10.049, 76.041
20.024, 90.146
50.272, 99.334
99.056, 99.922
NaN, NaN
];

d100 = [
1.0118, 49.786
1.9747, 53.022
4.9720, 55.360
10.079, 69.010
20.098, 81.379
50.014, 91.479
99.136, 98.004
200.23, 100.24
];

d200 = [
1.0053, 45.128
1.9950, 48.913
5.0109, 57.005
10.070, 71.110
20.067, 85.032
50.028, 90.840
99.098, 98.917
200.30, 99.416
];

override_x = [1,2,5,10,20,50,100,200];

z = [2, 5, 10, 50, 100, 200];

y = [d2(:,2), d5(:,2), d10(:,2), d50(:,2), d100(:,2), d200(:,2)];

x = override_x;

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    
    cmap = lcmp(numel(z));

    for i = 1:numel(z)
        l = z(i);
        
        plot(x, y(:,i), 's:', ...
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

