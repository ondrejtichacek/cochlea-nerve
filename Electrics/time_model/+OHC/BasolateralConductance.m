function [ G, dG, G_pre, dG_pre, dep, dep_pre, channels ] = BasolateralConductance( xgrid, simulation_units, args )
%BASOLATERALCONDUCTANCE
%OHC

arguments
    xgrid (:,1) double = linspace(0, 1, 300)'
    simulation_units (1,1) struct = struct()
    args.version = 'Johnson_2011'
    args.plotflag = false
    args.save_figures = false
end

%%

% [char_frequency, char_position] = frequency_map();
[char_position, char_frequency] = characteristic_frequency_model('Li');

default_tikz_options = { ...
    'relativeDataPath', 'img/HairCells/', ...
    'extraaxisoptions', {'legend style={font=\tiny}'}, ...
    'parseStrings', true};

sim_to_mV = Unit.conversionConstant(simulation_units.voltage, 'mV');
mV_to_sim = Unit.conversionConstant('mV', simulation_units.voltage);

%%

switch args.version
    
    case 'Tichacek_2022_v1'

        channels = BasolateralChannelsOHC;

        channels = channels([1,3]);

        % cf = [350, 2500, 12000]; % [Hz]
        % cx = char_position(cf);

        % distribute total conductance
        G1 = @(x) (25.8364 + boltzmann(x, 424.0532, 0.1592, -0.1513)) .* linf(0,1,x); % LF channels
        G2 = @(x) (25.8364 + boltzmann(x, 424.0532, 0.1592, -0.1513)) .* linf(1,0,x); % HF channels
        
        G1_pre = G1(xgrid);
        G2_pre = G2(xgrid);
        
        conductance = {@(x,P) G1(x) .* P, @(x,P) G2(x) .* P}; % nS
        conductance_der = {@(x) G1(x), @(x) G2(x)}; % nS/dP
        
        conductance_pre = {@(P) G1_pre .* P, @(P) G2_pre .* P}; % nS

        % TODO -- these functions are in fact constatnt - include in jac_0
        conductance_der_pre = {@(P) G1_pre, @(P) G2_pre}; % nS/dP
        
        G = conductance;
        dG = conductance_der;

        G_pre = conductance_pre;
        dG_pre = conductance_der_pre;
    
        dep = {{'xpos', 'popen'}, {'xpos', 'popen'}};
        dep_pre = {{'popen'}, {'popen'}};
        
        if args.plotflag == true
            hfig = figure;
            hold on
            xx = linspace(0, 1, 10);
            V = linspace(-120e-3, 0e-3, numel(xgrid));
            colors = ["b", "r"];
            for j = 1:numel(channels)
                for i = 1:numel(xx)
                    [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(j));            
                    plot(V*1e3, conductance{j}(xx(i), Oinf), ...
                        'Color', colors(j), ...
                        'DisplayName', channels(j).name);
                end
            end
            title('OHC basolateral conductance');
            xlabel('Voltage (mV)')
            ylabel('Conductance (nS)')
            legend('Location', 'NorthWest');
            
        end

    
    case 'dev_v1'
        
        channels = BasolateralChannelsOHC();
        
        % ~ half of conductance through slow, half through fast channels
        G1 = @(x) (25.8364 + boltzmann(x, 424.0532, 0.1592, -0.1513)) / 2;
        G2 = @(x) (25.8364 + boltzmann(x, 424.0532, 0.1592, -0.1513)) / 2;
        
        G1_pre = G1(xgrid);
        G2_pre = G2(xgrid);
        
        conductance = {@(x,P) G1(x) .* P, @(x,P) G2(x) .* P}; % nS
        conductance_der = {@(x) G1(x), @(x) G2(x)}; % nS/dP
        
        conductance_pre = {@(P) G1_pre .* P, @(P) G2_pre .* P}; % nS

        % TODO -- these functions are in fact constatnt - include in jac_0
        conductance_der_pre = {@(P) G1_pre, @(P) G2_pre}; % nS/dP
        
        G = conductance;
        dG = conductance_der;

        G_pre = conductance_pre;
        dG_pre = conductance_der_pre;
    
        dep = {{'xpos', 'popen'}, {'xpos', 'popen'}};
        dep_pre = {{'popen'}, {'popen'}};
        
        if args.plotflag == true
            hfig = figure;
            hold on
            xx = linspace(0, 1, 10);
            V = linspace(-120e-3, 0e-3, numel(xgrid));
            colors = ["b", "r"];
            for j = 1:numel(channels)
                for i = 1:numel(xx)
                    [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(j));            
                    plot(V*1e3, conductance{j}(xx(i), Oinf), ...
                        'Color', colors(j), ...
                        'DisplayName', channels(j).name);
                end
            end
            title('OHC basolateral conductance');
            xlabel('Voltage (mV)')
            ylabel('Conductance (nS)')
            legend('Location', 'NorthWest');
            
        end
    
    case 'LopezPoveda_2006'
    % Not OHC - tmp substitution
        
        channels = BasolateralChannelsOHC();
        
        G1 = channels(1).G * 1e9; % in nS;
        G2 = channels(2).G * 1e9; % in nS;
        
        conductance = {@(P) G1 * P, @(P) G2 * P}; % nS
        conductance_der = {G1, G2}; % nS/dP
        
        G = conductance;
        dG = conductance_der;

        G_pre = conductance;
        dG_pre = conductance_der;
    
        dep = {{'popen'}, {'popen'}};
        dep_pre = {{'popen'}, {'popen'}};
        
        if args.plotflag == true
            
            hfig = figure(3);
            hold on
            V = linspace(-120e-3, 0e-3, 200);
            for j = 1:numel(channels)
                [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(j));            
                plot(V*1e3, Oinf .* channels(j).G * 1e9, ...
                    'DisplayName', channels(j).name);
            end
            title('IHC basolateral conductance');
            xlabel('Voltage (mV)')
            ylabel('Conductance (nS)')
            legend('Location', 'NorthWest');
            
        end
        
    case 'Johnson_2011'
    % Data from Johnson 2011, Figure 5.E
    % Parameters for the Boltzmann function:
    % Voltage dependence of the K+ channel conductance
    
        channels = [];
        
        CF = [350, 900, 2500, 12000]; % [Hz]

        M = [30, 60, 107, 280]; % [nS]
        H = [-52, -63, -67, -56]; % [mV]
        S = [9.2, 8.4, 11, 18]; % [mV]

        if args.plotflag

            limits = [-120,0];
            V = linspace(limits(1), limits(2), 100);

            hfig = figure;
            hold on

            legend_text = cell(1,numel(CF));
            for i = 1:length(CF)
                plot(V, boltzmann(V, M(i), H(i), S(i)))
                legend_text{i} = sprintf('%g kHz', CF(i)/1000);
            end
            legend(legend_text, 'Location', 'NorthWest')

            title('OHC basolateral conductance');
            xlabel('V [mv]');
            ylabel('OHC G_K [nS]');

            if args.save_figures
                mySaveFig(hfig, 'OHC_G_K_fit', [], 'Results/HairCells/', default_tikz_options{:});
            end
        end
        
        extra_X = [  1,     0];
        extra_M = [ 30,   340];
        extra_H = [-48,   -52];
        extra_S = [  9,    20];

        % initial estimate for fitting
        x0M = [ 34.0848  259.1947    0.3025   -0.0838];
        x0H = [-46.4325   12.3787    0.3820   -0.1461];
        x0S = [ 19.4219  -11.0177    0.5684    0.0864];

        X = char_position(CF);

        [FM, xbest, fval] = fitting.boltzmann([X, extra_X], [M, extra_M], x0M);
        [FH, xbest, fval] = fitting.boltzmann_der([X, extra_X], [H, extra_H], x0H);
        [FS, xbest, fval] = fitting.boltzmann([X, extra_X], [S, extra_S], x0S);

        % x0H = [-67  14  0.25 -0.0838  15  0.6 0.04];
        % FUN = @(t,x1,x2,x3,y1,y2,y3) boltzmann(t,x1,x2,x3) + boltzmann(t,y1,y2,y3);
        % do_fit = false;
        % [FH, xbest, fval] = fitting.optim_fit(FUN, [X, extra_X], [H, extra_H], x0H, [], [], do_fit);

        if args.plotflag
            xx = linspace(0, 1, 100);

            hfig = figure;

            subplot(1,3,1)
            hold on
            plot(X, M, 'LineStyle', 'none', 'marker', 'o', 'Color', 'b');
            plot(extra_X, extra_M, 'LineStyle', 'none', 'marker', 'x', 'Color', 'r');
            plot(xx, FM(xx), 'LineStyle', '--', 'Color', 'k');
            xlabel('Distance from stapes (norm.)');
            ylabel('M [nS]');

            subplot(1,3,2)
            hold on
            plot(X, H, 'LineStyle', 'none', 'marker', 'o', 'Color', 'b');
            plot(extra_X, extra_H, 'LineStyle', 'none', 'marker', 'x', 'Color', 'r');
            plot(xx, FH(xx), 'LineStyle', '--', 'Color', 'k');
            xlabel('Distance from stapes (norm.)');
            ylabel('H [mV]');

            subplot(1,3,3)
            hold on
            plot(X, S, 'LineStyle', 'none', 'marker', 'o', 'Color', 'b');
            plot(extra_X, extra_S, 'LineStyle', 'none', 'marker', 'x', 'Color', 'r');
            plot(xx, FS(xx), 'LineStyle', '--', 'Color', 'k');
            xlabel('Distance from stapes (norm.)');
            ylabel('S [mV]');

            %
            % plot new curves from interpolated data
            xxx = linspace(0, 1, 10);

            if args.save_figures
                mySaveFig(hfig, 'OHC_G_K_fit_parameters', [], 'Results/HairCells/', default_tikz_options{:}, ...
                    'width', '\textwidth');
            end

            hfig = figure;
            hold on

            cmp = rwb(numel(CF));
            legend_text = cell(1,numel(CF));
            for i = 1:length(CF)
                plot(V, boltzmann(V, M(i), H(i), S(i)), 'Color', cmp(i,:), 'LineStyle', '-', 'Marker', 'none')
                legend_text{i} = sprintf('%g kHz', CF(i)/1000);
            end

            cmp = flip(rwb(numel(xxx)));
            legend_text2 = cell(1,numel(xxx));
            for i = 1:length(xxx)
                plot(V, boltzmann(V, FM(xxx(i)), FH(xxx(i)), FS(xxx(i))), 'Color', cmp(i,:), 'LineStyle', '-.', 'Marker', 'none')
                legend_text2{i} = sprintf('%g kHz', char_frequency(xxx(i))/1000);
            end

        %     legend([legend_text, legend_text2], 'Location', 'NorthWest')
            legend(legend_text, 'Location', 'NorthWest')

            % set(gca,'YScale', 'log')

            title('OHC basolateral conductance');
            xlabel('V [mv]');
            ylabel('OHC G_K [nS]');

            if false
                cbarTL = cellfun(@(x) sprintf('%d', round(Frequency(char_frequency(x), 'Hz').Hz)), num2cell(xxx), 'UniformOutput', false);

                manual_colorbar(cmp, cbarTL, ...
                    'orientation', 'vertical', ...
                    'reversed', true, ...
                    'label', sprintf('f (Hz)'));
            end

            if args.save_figures
                mySaveFig(hfig, 'OHC_G_K_fit_extrapolation', [], 'Results/HairCells/', default_tikz_options{:});
            end
        end
        
        G = @(x,V) boltzmann(V * sim_to_mV, FM(x), FH(x), FS(x));
        dG = @(x,V) boltzmann_der(V * sim_to_mV, FM(x), FH(x), FS(x)) / mV_to_sim;

        M = FM(xgrid); H = FH(xgrid); S = FS(xgrid);

        G_pre = @(V) boltzmann(V * sim_to_mV, M, H, S);
        dG_pre = @(V) boltzmann_der(V * sim_to_mV, M, H, S) / mV_to_sim;
        
        dep = [];
        dep_pre = [];
        
    otherwise
        error('Unknown version %s', args.version)
end


    function [b, a] = linear_regression(X, y)
        
        beta = X\y;
        
        a = beta(2);
        b = beta(1);
    end

end

