function [ hfig ] = profile_one_graph( xx, volt_fun, variables, mnaopt, plotopt, runopt, analysisStartFrame )
%PROFILE_ONE_GRAPH 

for relative = {true, false}
    hfig = plotopt.figure;
    % hfig = figure(100);

    hold on

    legendtext = cell(1, numel(variables));
    c = rwb(numel(variables));

    for i = 1:numel(variables)
        H(i) = plot_steady_state_and_profile(xx, ...
            volt_fun(i), analysisStartFrame, 1, relative{:}, 'lineProps', {'Color', c(i,:)});
        legendtext{i} = variables{i};
    end
    l = legend([H.mainLine], legendtext, 'Location', 'SouthEast');
    l.ItemHitFcn = @plotOpt.emph_line_width;

    xlabel('BM position');
end
end

