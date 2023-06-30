function [x, y, z] = fig9a(args)
%FIG9a
arguments
    args.plotflag (1,1) logical = false
end

p10 = [
    0.99445, 86.074
1.9895, 82.905
4.9086, 80.873
9.9602, 85.711
19.857, 89.025
49.718, 99.575
NaN, NaN
NaN, NaN
];

p20 = [
0.99477, 55.188
1.9917, 58.024
4.8951, 59.805
9.9794, 67.122
19.917, 78.920
49.491, 97.955
98.221, 100.03
294.86, 96.656

];

p30 = [
0.98178, 21.442
1.9911, 22.371
4.8962, 28.252
9.9425, 38.809
NaN, NaN
49.708, 64.590
97.513, 78.104
293.22, 87.219

];

p40 = [
0.99350, 12.099
1.9973, 12.457
4.8887, 16.526
9.9226, 23.462
NaN, NaN
49.601, 48.193
97.339, 64.472
292.88, 78.353

];

p60 = [
0.98844, 6.3794
NaN, NaN
4.9032, 5.7539
9.9513, 12.213
NaN, NaN
49.548, 39.995
97.248, 57.322
292.69, 73.491

];



override_x = [1,2,5,10,20,50,100,300];

z = [10, 20, 30, 40, 60];

y = [p10(:,2), p20(:,2), p30(:,2), p40(:,2), p60(:,2)];

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

