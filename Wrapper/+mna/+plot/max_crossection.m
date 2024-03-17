function [ ] = max_crossection( time, T, t0, tf, IHC_V, mnaopt, plotopt, runopt )
%MAX_CROSSECTION 

if plotopt.subplot == true
    f3 = plotopt.figure;
else
    sf9 = plotopt.figure;
end


plot(time.ms, IHC_V, 'Color', [.5 .5 .5]);
hold on

IHC_V = interp1(time.ms, IHC_V, T);

plot(T, IHC_V, '-bo');

WindowLength = round(Time(1,'ms') * mnaopt.samplingFrequency);
Y = movingAverageFilter(WindowLength, IHC_V);
plot(T, Y, 'Color', 'k', 'LineWidth', 0.75);

title('Maximal IHC crossection');
xlabel('Time, [ms]');
ylabel('IHC Voltage, [mV]');
xlim([t0.ms, tf.ms]);

%     h = plotopt.figure;
%     plot(TIME*1e3, IHC_V);
%     title('Maximal IHC crossection');
%     xlabel('Time, [ms]');
%     xlim([0 TIME(end)*1e3]);
%     ylabel('IHC Voltage, [mV]');
%     
% save figures
if runopt.save_figures
    if plotopt.subplot == true
        plotopt.save_figure(f3,  runopt.figureSaveDir, 'figure_3');
%             plotopt.save_figure(h,  runopt.figureSaveDir, 'MaxCrossectionIHC'));

    else
        plotopt.save_figure(sf9,  runopt.figureSaveDir, 'sfigure_9');
    end
end


end

