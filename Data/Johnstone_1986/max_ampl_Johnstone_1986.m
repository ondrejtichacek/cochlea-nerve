function [] = max_ampl_Johnstone_1986(plotflag, runopt)
arguments
    plotflag (1,1) logical = false
    runopt = []
end
%MAX_AMPL_JOHNSTONE_1986 Summary of this function goes here
%   Detailed explanation goes here

% Johnstone et al Hear Res. 1986;22:147-53. doi: 10.1016/0378-5955(86)90090-0

% maximal amplitude vs dB SPL (normal)
% column 1 ... db SPL
% column 2 ... amplitude [nm]
max_ampl_Johnstone_1986 = [ ...
    9.8970, 0.82694
    19.167, 3.1623
    29.373, 8.2233
    39.745, 13.373
    49.730, 20.335
    59.310, 21.029
    80.198, 23.648
    90.002, 43.975
    99.790, 60.472
    110.20, 202.22
]; % data from the original paper processed through http://arohatgi.info/WebPlotDigitizer/app/?

% maximal amplitude vs dB SPL (after trauma)
% column 1 ... db SPL
% column 2 ... amplitude [nm]
max_ampl_Johnstone_1986_trauma = [ ...
    50.799, 0.84800
    59.822, 6.7247
    70.188, 9.8888
    82.679, 19.997
    90.268, 29.904
    120.80, 586.41
]; % data from the original paper processed through http://arohatgi.info/WebPlotDigitizer/app/?


if plotflag
    hfig = figure();
    hold on

    % Experimental data -- normal
    plot(max_ampl_Johnstone_1986(:,1), max_ampl_Johnstone_1986(:,2), ...
        'LineStyle', 'none', ...
        'Marker', 's', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', 'k', ...
        'Color', 'k');

    % Experimental data -- trauma
    plot(max_ampl_Johnstone_1986_trauma(:,1), max_ampl_Johnstone_1986_trauma(:,2), ...
        'LineStyle', 'none', ...
        'Marker', 'o', ...
        'MarkerSize', 6, ...
        ...'MarkerFaceColor', 'k', ...
        'Color', 'k');
    
    plot([50,121], [0.9, 586.41], ...
        'Color', 'red');

    plot([5,40], [0.9, 20], ...
        'Color', 'blue');
    
    plot([30,100], [10, 50], ...
        'Color', 'black');
    
    xlim([0,125]);
    
    xlabel('Level [dB SPL]');
    ylabel('BM displacement [nm]');
    
    ax = gca();
    ax.YScale = 'log';
    
    legend('normal', 'trauma', 'location', 'SE');

    %% Save figure

    FIGNAME = 'max_ampl_Johnstone_1986';
    PATH = 'img/MECH/';

    if ~isempty(runopt) && runopt.save_figures
        analysis.plot.export_figure(hfig, FIGNAME, PATH, false, runopt)
    end
end

end

