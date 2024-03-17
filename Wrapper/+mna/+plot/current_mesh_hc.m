function [ ] = current_mesh_hc( Xgrid, Tgrid, CurrIHC, CurrOHC, t0, tf, plotopt, runopt, unit )
%CURRENT_MESH_HC 

if plotopt.subplot == true
    f1 = plotopt.figure;
end

% ---------------------------------------------------------------------
% IHC
if plotopt.subplot == true
    subplot(1,2,1);
else
    sf1 = plotopt.figure;
end

surf(Xgrid, Tgrid, CurrIHC, ... % from mA to pA
    plotopt.surf_options{:});

xlabel('BM position');
ylim([t0.ms, tf.ms]);
ylabel('Time [ms]');
zlabel(sprintf('Current (%s)', unit));
title('single IHC');
colormap(plotopt.colormap)
view(plotopt.view)

% ---------------------------------------------------------------------
% OHC
if plotopt.subplot == true
    subplot(1,2,2)
else
    sf2 = plotopt.figure;
end


surf(Xgrid, Tgrid, CurrOHC, ...
    plotopt.surf_options{:});

xlabel('BM position');
ylim([t0.ms, tf.ms]);
ylabel('Time [ms]');
zlabel(sprintf('Current (%s)', unit));
title('single OHC');
colormap(plotopt.colormap)
view(plotopt.view)

% ---------------------------------------------------------------------
% save figures
if runopt.save_figures
    if plotopt.subplot == true
        plotopt.save_figure(f1,  runopt.figureSaveDir, 'figure_1');
    else
        plotopt.save_figure(sf1,  runopt.figureSaveDir, 'sfigure_1');
        plotopt.save_figure(sf2,  runopt.figureSaveDir, 'sfigure_2');
    end
end

end

