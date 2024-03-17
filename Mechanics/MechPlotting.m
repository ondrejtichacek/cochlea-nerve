function [status] = MechPlotting(resFiles, plotopt, stimulus, topt, midopt, mechopt, runopt, opt, memopt, paropt)
%MECHPLOTTING
arguments
    resFiles (1,1) struct
    plotopt (1,1) plotOpt
    stimulus (1,1) {mustBeA(stimulus, ["Signal", "CompoundSignal"])}
    topt (1,1) timeOpt
    midopt (1,1) midOpt
    mechopt (1,1) mechOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt = parOpt()
end

status = true;

if ~plotopt.doplot
    return
end


doplot = plotopt.do.mechanical;

if islogical(doplot) && ~doplot
    return
end

if ~any(structfun(@(x) x, doplot))
    return
end

%% clear stored figure handles

plotopt.figure_handles = [];

%%


fBM   = resFiles.fBM;
% fmech = resFiles.fmech;

x = fBM.x;

%% undersampling

fsampldesired = Frequency(20, 'kHz');

%% TIME

time = Time(fBM.time, fBM.time_unit);

[T, ts, t0, tf] = plotopt.time_samples(time, topt, fsampldesired);

[Xgrid,Tgrid] = meshgrid(x, T);


%% plots

variables = {'BMx', 'BMv', 'TMx', 'TMv', 'TMf'};
descriptions = {'BM amplitude', 'BM velocity', 'stereocilia deflection', 'stereocilia deflection velocity', 'stereocilia filtered'};
units = {'nm', '', '', '', ''};

for i = 1:numel(variables)

    if doplot.(variables{i})

        hfig = plotopt.figure;

        VAR = plotopt.subsample(fBM, variables{i}, ts, ':');

        % VAR = BM_saturation(VAR);
        
        M = max(abs(VAR(:)));
        
        surf(Xgrid, Tgrid, VAR, ...
            plotopt.surf_options{:});

        caxis([-M M]);
        
    %     hold on
    %     mesh(x, TIME, VAR, ...
    %         'FaceColor','interp', ...
    %         ...'EdgeColor','k', 'EdgeAlpha', 0.1 ...
    %         'EdgeColor','none' ...
    %         );
    %     mesh(x, TIME2, VAR, ...
    %         'FaceColor','none', ...
    %         'EdgeColor','k', 'EdgeAlpha', 0.13, 'MeshStyle', 'row' ...
    %         ...'EdgeColor','none' ...
    %         );

        xlabel('BM position');
        ylim([t0.ms, tf.ms]);
        ylabel('Time [ms]');
        zlabel(sprintf('%s %s',descriptions{i}, units{i}));
        colormap(plotopt.centered_colormap)
        view(plotopt.view)

        if runopt.save_figures
            plotopt.save_figure(hfig,  runopt.figureSaveDir, variables{i});
        end

    end

    if doplot.([variables{i},'_profile'])

        hfig = plotopt.figure;

        plot_steady_state_and_profile(x, fBM.(variables{i}));

        xlabel('BM position');
        ylabel(sprintf('%s %s',descriptions{i}, units{i}));
        title(sprintf('%s profile',descriptions{i}));

        if runopt.save_figures
            %plotopt.save_figure(hfig,  runopt.figureSaveDir, sprintf('ss_%s', variables{i}));
            mySaveFig(hfig, sprintf('ss_%s', variables{i}), runopt.path.mech, runopt.path.mech, ...
                'height', '0.2\textwidth', ...
                'width', '0.4\textwidth', ...
                'showInfo', false);
        end

    end

end

%% Show figures

if ~isempty(plotopt.figure_handles)
    arrayfun(@(h) set(h, 'Visible', 'On'), plotopt.figure_handles);
end


end
