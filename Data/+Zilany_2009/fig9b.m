function [x, y, z] = fig9b(args)
%FIG9b
arguments
    args.plotflag (1,1) logical = false
end

data_onset = [
1.7764e-15, 205.17
5.0239, 188.98
14.990, 204.20
20.014, 161.85
25.017, 230.15
30.021, 235.02
35.024, 215.32
];

data_window = [
0.020506, 135.32
5.0034, 149.56
15.010, 154.63
20.014, 155.41
24.997, 158.15
29.979, 148.00
35.003, 140.78
];

model_onset = [
-0.14354, 251.41
4.9624, 228.98
15.113, 239.32
20.014, 198.93
25.058, 212.59
30.103, 226.05
35.106, 220.20
];

model_window = [
0.0000, 103.90
5.0239, 113.85
15.010, 114.24
20.014, 114.05
24.997, 106.63
30.000, 137.27
35.003, 120.29
];



override_x = [0,5,15,20,25,30,35];

y = [data_onset(:,2), data_window(:,2), model_onset(:,2), model_window(:,2)];

x = override_x;

if args.plotflag
    figure
    hold on
    set(gca, 'xscale', 'log')
    
   
    plot(x, y(:,1), 'o-', ...
        'Color', 'b', ...
        'DisplayName', 'onset, experiment, Smith et al. (1985)')
   
    plot(x, y(:,2), 'x-', ...
        'Color', 'b', ...
        'DisplayName', 'window, experiment, Smith et al. (1985)')

    plot(x, y(:,3), 'o:', ...
        'Color', [0.5,0.5,0.5], ...
        'DisplayName', 'onset, model, Zilany et al. (2009)')
   
    plot(x, y(:,4), 'x:', ...
        'Color',  [0.5,0.5,0.5], ...
        'DisplayName', 'window, model, Zilany et al. (2009)')    

    xlabel('Time Delay (ms)')
    ylabel('Increment in spikes/s')
    legend('Location', 'best')

    %xlim([0.9,330])
    %ylim([10, 105])
end

end

