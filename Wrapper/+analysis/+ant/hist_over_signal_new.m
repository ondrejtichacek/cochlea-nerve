function [data] = hist_over_signal_new(RES, SIGNAL_TO_PLOT, ...
        Configurations, Conf_nums, Conf_info,...
        runopt, opt, num_replicas, stimulus_fcn, args)
%HIST_OVER_SIGNAL
arguments
    RES
    SIGNAL_TO_PLOT
    Configurations
    Conf_nums
    Conf_info
    runopt
    opt
    num_replicas
    stimulus_fcn
    args.BinWidth (1,1) Time = Time(2, 'ms')
    args.plotflag (1,1) logical = true
    args.hfig = []
    args.conf_DisplayName = []
end


runopt.update('figureSaveDir', fullfile(opt.cochleadir, 'Results', 'PP'));

data = cell(Conf_nums);

if args.plotflag
    if isempty(args.hfig)
        hfig = figure();
    else
        hfig = args.hfig;
    end
end

cmap = lcmp(numel(Configurations));

for i = 1:numel(Configurations)
    unit = 's';
    Location{i} = Time(RES.(SIGNAL_TO_PLOT).Location{i}, unit);

    t = args.BinWidth;
    t2 = 3*t;

    [counts{i}, edges{i}] = histcounts( Location{i}.ms, ...
        'BinWidth', t.ms, ...
        'Normalization', 'count');

    [counts_2{i}, edges_2{i}] = histcounts( Location{i}.ms, ...
        'BinWidth', t2.ms, ...
        'Normalization', 'count');

    ff = 1/t;
    f2 = 1/t2;

    n_dyn_rep =  100;

    counts{i} = ff.Hz * counts{i} / num_replicas / n_dyn_rep;
    counts_2{i} = f2.Hz * counts_2{i} / num_replicas / n_dyn_rep;
end

max_counts = max(cat(2, counts{:}));

for i = 1:numel(Configurations)

    if args.plotflag
        figure(hfig)
        nexttile
        hold on

        histogram( ...
            'BinEdges', edges{i}, ...
            'BinCounts', counts{i}, ...
            'EdgeColor', cmap(i,:), ...
            'FaceColor', cmap(i,:), ...
            ...'FaceColor', 'k', ...
            'FaceColor', cmap(i,:), ...
            'FaceAlpha', 1)

        stairs(edges_2{i}(1:end-1), counts_2{i}, 'k-')
    
        if ~isempty(args.conf_DisplayName)
            text(0, max_counts*0.7, args.conf_DisplayName{i})
        end

        % histogram( Location.ms, ...
        %     'BinWidth', t.ms, ...
        %     'FaceColor', 'k', ...
        %     'FaceAlpha', 1, ...
        %     'Normalization', 'count' );
        
        %  ylim([0,180])
        xlabel('Time (ms)');
        ylabel('');

        ylim([0, ceil(max_counts)])
    end

    Conf = Configurations(i);
    [~, stimulus] = stimulus_fcn(Conf);

    M = max(counts{i});
    % MM(r,s,j,i) = M;

    if args.plotflag && ismethod(stimulus, 'envelope_fcn')
        envelope = stimulus.envelope_fcn('ms');

        plot(stimulus.audiotime.ms, M * envelope(stimulus.audiotime.ms), 'r:')
    end

    % if runopt.save_figures
    %     fname = sprintf('Peaks_Time_hist_%s_%s_%d', SR_TO_PLOT{s}, SIGNAL_TO_PLOT{r}, FREQUENCY(j));
    %     mySaveFig(hfig, fname, [], runopt, 'width', '0.9\textwidth', 'height', '0.2\textwidth', 'showInfo', false);
    % end

    t = Time(0.5, 'ms');

    if isempty(Location{i}.ms)
        xi = [];
        f = [];
    else
        [f, xi] = ksdensity(Location{i}.ms, 'Bandwidth', t.ms);
        
        if args.plotflag
            % plot(xi, M * f / max(f), 'k')
        end
    end

    data{i} = struct( ...
        counts=counts{i}, ...
        edges=edges{i}, ...
        kds_x=xi, ...
        kds_f=f);

end

end

