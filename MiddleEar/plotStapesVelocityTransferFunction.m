function [hfig] = plotStapesVelocityTransferFunction(FREQUENCY, Y, P, args)
%PLOTTRANSFERFUNCTION
arguments
    FREQUENCY = []
    Y = []
    P = []
    args.hfig (1,:) = [figure, figure];
    
    args.do_plot_all (1,1) logical = true
    args.do_plot_selected (1,:) cell = {}
    
    args.plot_magnitude (1,1) logical = true %false
end

if not(isempty(args.do_plot_selected))
    args.do_plot_all = false;
end

in = @(s, c) any(strcmp(s, c));

hfig = args.hfig;

cmap = parula(16);

col = struct( ...
    'green', cmap(10,:), ...
    'yellow', cmap(14,:));

%% Velocity

figure(hfig(1));
if numel(hfig) == 1
    subplot(2,1,1)
end
hold on

if ~args.plot_magnitude
    set(gca, 'YScale', 'log')
end

if args.do_plot_all || in('Chien_2009', args.do_plot_selected)
    
    if args.plot_magnitude
        [x, ~, ~, y, yerr] = Chien_2009.fig9a();        
    else
        [x, y, yerr] = Chien_2009.fig9a();
    end
    
    shadedErrorBar(x, y, yerr, ...
        'lineProps', {'Color', col.yellow, 'DisplayName','Chien (2009)'});
    
    % if args.plot_magnitude
    %     [x, ~, ~, y, yerr] = Chien_2009.fig10a_human_live();
    % else
    %     [x, y, yerr] = Chien_2009.fig10a();
    % end
    % 
    % shadedErrorBar(x, y, yerr, ...
    %     'lineProps', {'Color', col.yellow, 'DisplayName','Chien (2009)'});
    % 
    % if args.plot_magnitude
    %     [x, ~, y] = Chien_2009.fig10a_human_cadaveric();
    % else
    %     [x, y, ~] = Chien_2009.fig10a_human_cadaveric();
    % end
    % 
    % plot(x, y, 'Color', col.yellow, 'DisplayName','Chien (2009)')
    
end
    
if args.do_plot_all || in('Aibara_2001', args.do_plot_selected)
    
    if args.plot_magnitude
        [x, ~, ~, y, yerr] = Aibara_2001.fig7a();
    else
        [x, y, yerr] = Aibara_2001.fig7a();
        
        % Aibara gives the velocity in 'mm/s/Pa' - convert to um/s/Pa
        k = 1e3;

        % This is needed shomehow. Probably different quantities?
        k =  k * 0.18;

        y = y * k;
        yerr = yerr * k;
    end
    
    
    shadedErrorBar(x, y, yerr, ...
        'lineProps', {'Color', col.green, 'DisplayName','Aibara (2001)'});
end

if ~isempty(FREQUENCY)
    plot(FREQUENCY, Y, ...
        'Color', 'k', 'LineWidth', 2, 'DisplayName','model');
end

set(gca, 'XScale', 'log')
erbticks;
xlabel('Frequency (Hz)');
ylabel('Velocity (\mum/s/Pa)');
grid on
% xticks([125, 250, 500, 1000, 2000, 4000, 8000, 16000]);
legend('Location', 'South')

%% Phase lag

if numel(hfig) == 1
    figure(hfig(1));
    subplot(2,1,2)
else
    figure(hfig(2));
end
hold on

if args.do_plot_all || in('Chien_2009', args.do_plot_selected)
    [x, y, yerr] = Chien_2009.fig9b();
    
    % Chien gives the phase in 'periods' - convert to degrees
    k = 360;
    
    y = y * k;
    yerr = yerr * k;
    
    shadedErrorBar(x, y, yerr, ...
        'lineProps', {'Color', col.yellow, 'DisplayName','Chien (2009)'});
    
end

if args.do_plot_all || in('Aibara_2001', args.do_plot_selected)
    [x, y, yerr] = Aibara_2001.fig7b();
    shadedErrorBar(x, y, yerr, ...
        'lineProps', {'Color', col.green, 'DisplayName','Aibara (2001)'});
end

if ~isempty(FREQUENCY)
    plot(FREQUENCY, P, ...
        'Color', 'k', 'LineWidth', 2, 'DisplayName','model');
end

set(gca, 'XScale', 'log')
erbticks;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on
% xticks([125, 250, 500, 1000, 2000, 4000, 8000, 16000]);
% yticks([-270, -180, -90, 0, 90, 180]);
legend('Location', 'South')

end

