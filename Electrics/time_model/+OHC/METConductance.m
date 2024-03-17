function [ MET_conductance, MET_conductance_der, MET_conductance_pre, MET_conductance_der_pre ] = METConductance( xgrid, simulation_units, args )
%METCONDUCTANCE
%OHC

arguments
    xgrid (:,1) double = linspace(0, 1, 300)'
    simulation_units (1,1) struct = struct()
    args.version = 'Johnson_2011'
    args.plotflag = false
    args.save_figures = false

    args.noise_damage (1,1) logical = false
    args.noise_damage_f0 (1,1) Frequency = Frequency(1200, 'Hz')
    args.noise_damage_degree (1,1) double = 0.1
end

%%

% [char_frequency, char_position] = frequency_map();
[char_position, char_frequency] = characteristic_frequency_model('Li');

default_tikz_options = { ...
    'relativeDataPath', 'img/HairCells/', ...
    'extraaxisoptions', {'legend style={font=\tiny}'}, ...
    'parseStrings', true};

%%

switch args.version

    case 'Johnson_2011'
        
        % -----------------------------------------------------------------
        % MET open probability

        MET_params_cilia = {8*0.5/30, 35*0.5/30};
        MET_popen_cilia = @(y) boltzmann(y, 1, MET_params_cilia{:});
        
        MET_popen_cilia_mag = @(y, m) boltzmann(y, m, MET_params_cilia{:});
        MET_popen_cilia_mag_der = @(y, m) boltzmann_der(y, m, MET_params_cilia{:});
        
        if args.plotflag
    
            % TM = linspace(-80, 80, 500);
            TM = linspace(-6, 6, 500);
            TM0 = 0;

            hfig = figure;
            hold on

            plot(TM, MET_popen_cilia(TM), 'Color', 'k', 'LineWidth', 0.75);
            plot(TM0, MET_popen_cilia(TM0), 'Color', 'k', 'Marker', '*');
            plot([TM(1), TM0], [MET_popen_cilia(TM0), MET_popen_cilia(TM0)], ...
                'Color', [0.7, 0.7 0.7], ...
                'LineStyle', ':', ...
                'LineWidth', 0.75);

            xlabel('Stereocilia displacement [nm]')
            ylabel('p_{open}')
            title('OHC MET open probability')

            if args.save_figures
                mySaveFig(hfig, 'MET_OHC', [], 'Results/HairCells/', default_tikz_options{:});
            end

        end
        
        % -----------------------------------------------------------------
        % max MET conductance
        % Data from Johnson 2011, Supelementary Figure 1.B
        
        fig = Johnson_2011_data();
        
        [~, ~, cf, cond] = combineDatasets(fig(9).B.datasets, 'rat', 'gerbil');

        % OHC.MET_conductance_max = @(x) 20 + 70*exp(-x/0.3);
        FUN = @(t, x1, x2) x1.*exp(-t./x2);
        MET_conductance_max = fitting.optim_fit(FUN, char_position(cf(:)), cond(:), [20, 70, 0.3]);

        if args.plotflag
            xx = linspace(0,1,100);    

            hfig = figure;
            hold on

            legendtext = Johnson_2011_plot(fig(9), 'B', char_position);

            set(gca,'XScale', 'linear')
            xlabel('Distance from stapes (norm.)')

            plot(xx, MET_conductance_max(xx), 'k');

            legendtext{end+1} = 'fit';

            title('OHC MET conductance')
        end

        if false
        % [F, xbest, fval] = fitting.boltzmann(cf(:), log10(cond(:)), [0.9258    1.8985    4.2876    1.4181], [], @log10)

        % OHC.MET_conductance_max = @(cf) 20 + 10.^F(cf);

        MET_conductance_max = @(cf) 20 + boltzmann(log10(cf)/log10(20000), 220, 1.1, 0.13);

        if args.plotflag
            cf = logspace(log10(50), log10(20000), 100);

            hfig = figure;
            hold on

            legendtext = Johnson_2011_plot(fig(9), 'B');
            plot(cf, OHC.MET_conductance_max(cf), 'LineStyle', '--', 'Color', 'k');
            set(gca,'XScale', 'log')

            title('OHC G_{MT}');
            xlabel('CF (Hz)');
            ylabel('G [nS]');


            if args.save_figures
                mySaveFig(hfig, 'OHC_MET_conductance', [], 'Results/HairCells/', default_tikz_options{:});
            end
        end
        end
        
        G_max_pre = MET_conductance_max(xgrid);

        if args.noise_damage
            [damage_fac, damage_fun] = noise_damage(args.noise_damage_f0, xgrid(:)', args.noise_damage_degree);
                
            G_max_pre = G_max_pre .* damage_fac(:);

            MET_conductance = @(x,y) MET_popen_cilia_mag(y, MET_conductance_max(x) .* damage_fun(x));
            MET_conductance_der = @(x,y) MET_popen_cilia_mag_der(y, MET_conductance_max(x) .* damage_fun(x));
        else
            MET_conductance = @(x,y) MET_popen_cilia_mag(y, MET_conductance_max(x));
            MET_conductance_der = @(x,y) MET_popen_cilia_mag_der(y, MET_conductance_max(x));
        end
                
        MET_conductance_pre = @(y) MET_popen_cilia_mag(y, G_max_pre);
        MET_conductance_der_pre = @(y) MET_popen_cilia_mag_der(y, G_max_pre);
        
        
    otherwise
        error('Unknown version %s', args.version)
end

    function [b, a] = linear_regression(X, y)
        
        beta = X\y;
        
        a = beta(2);
        b = beta(1);
    end

end

