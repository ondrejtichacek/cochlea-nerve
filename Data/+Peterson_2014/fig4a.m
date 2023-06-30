function [x, y] = fig4a(args)
%FIG4a
arguments
    args.plotflag (1,1) logical = false
end

xy = [ ...
    1.0000, 0.95398
    1.9945, 0.94281
    4.0370, 0.91545
    8.0520, 0.87848
    15.826, 0.85299
    31.104, 0.81855
    62.039, 0.78088
    123.74, 0.72333
    307.64, 0.78088
    613.59, 0.89414
    1241.9, 1.2145
    2477.1, 1.6496
    5013.7, 2.9379
    10000, 4.5963
];



x = xy(:,1);
y = xy(:,2);

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    
    plot(x, y(), 'ko-', ...
        'DisplayName', 'Peterson 2014')
   
    xlabel('Counting Time T (ms)')
    ylabel('Fano Factor F(T)')
    legend('Location', 'best')

    %xlim([0.9,330])
    %ylim([10, 105])
end

end

