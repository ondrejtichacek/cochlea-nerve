function [ ] = voltage_mesh( Xgrid, Tgrid, variables, ...
    t0, tf, plotopt, runopt )
%VOLTAGE_MESH 

if plotopt.subplot == true
    hfig = plotopt.figure;
end


if plotopt.subplot == true
    tiledlayout('flow');
else
    hfig = [];
end


names = fieldnames(variables);
for i = 1:numel(names)

    if plotopt.subplot == true
        nexttile
    else
        hfig(i) = plotopt.figure;
    end
    
    V = variables.(names{i});
    C = V.data - V.steadystate;
    
    M = max(abs(C(:)));
    
    surf(Xgrid, Tgrid, V.data, C, ...
            plotopt.surf_options{:});

    xlabel('BM position');
    ylabel('Time [ms]');
    zlabel('Voltage [mV]');
    
    ylim([t0.ms, tf.ms]);
    
    title(V.name);
    
    colormap(plotopt.centered_colormap)    
    caxis([-M, M])
    
    view(plotopt.view);
    
    set(gca, 'color', 0.95*ones(1,3))
end

% save figures
if runopt.save_figures
    if plotopt.subplot == true
        plotopt.save_figure(hfig, runopt.figureSaveDir, 'figure_2');
    else
        for i = 1:numel(hfig)
            plotopt.save_figure( ...
                hfig(i), runopt.figureSaveDir, sprintf('sfigure_%d', i));
        end
    end
end

end
