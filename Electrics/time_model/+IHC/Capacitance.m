function [ capacitance_apical, capacitance_basal, hfig ] = Capacitance( xgrid, simulation_units, args )
%CAPACITANCE
%IHC

arguments
    xgrid (:,1) double = linspace(0, 1, 300)'
    simulation_units (1,1) struct = struct()
    args.version = 'Johnson_2015'
    args.plotflag = false
    args.save_figures = false
end

%%

[char_frequency, char_position] = frequency_map();

default_tikz_options = { ...
    'relativeDataPath', 'img/HairCells/', ...
    'extraaxisoptions', {'legend style={font=\tiny}'}, ...
    'parseStrings', true};

%%

IHC = struct();

switch args.version

    case 'Johnson_2015'
        
        % Data from Johnson et al 2011 and Johnson 2015
        
        fig = Johnson_2011_data();

        % Johnson 2015 measured basal (~30 kHz) and apical (~300 Hz) turn in gerbil
        % IHCs.
        % hf ... high-freqency ... apical part of the cochlea
        % lf ... low-freqency ... basal part of the cochlea
        
        x_hf = 0; % 30 kHz is outside range in human
        x_lf = char_position(300);
        
        % we assume linear dependence on position
        % and perform linear regression
        % y = b0 + b1 * x
        xx = [x_hf; x_lf];
        X = [ones(length(xx),1), xx];
        
        C_hf = 11.7; % +- 0.2 pF % Johnson eLife 2015 (pp. 17-18)
        C_lf = 12.4; % +- 0.1 pF % Johnson eLife 2015 (pp. 17-18)

        IHC.Johnson_2015.capacitance = struct( ...
            'x', [x_hf, x_lf], ...
            'C', [C_hf, C_lf], ...
            'yerr', [0.2, 0.1]);

        [b, a] = linear_regression(X, [C_hf; C_lf]);
        
        IHC.capacitance_fit.linear = @(x) a*x + b;
        
    case 'Johnson_2011'
        
        
        IHC.Johnson_2011.capacitance = 12.5; % +- 0.5 pF % Johnson 2011 (pp. 1149)
        
        
    otherwise
        error('Unknown version %s', args.version)
end

%% Constant capacitance

% Johnson 2011 claims, that the ratio CA/CB of apical and basilar
% capacitance is r = 0.2 uniformly across the BM
% CA = CB*r, CB = CA/r
%
% For IHC we assume that the ratio is similar and we tune it to calibrate
% the time constants

capacitance_ratio = 0.1;

% Also, we know that the total membrane capacitance CM = CA + CB
% Therefore
% CA = CM / (1 + 1/r);
% CB = CM / (1 + r);

CA_factor = 1/(1 + 1/capacitance_ratio);
CB_factor = 1/(1 + capacitance_ratio);

% The OHC capacitance is not uniform along the tonotopy axis

% extract the function to prevent hashing error
% OHC_capacitance = @(x) IHC.capacitance_fit.linear(x);
capacitance_tmp = IHC.capacitance_fit.linear;
capacitance = @(x) capacitance_tmp(x);

% CA_pre = CA_factor * capacitance(xgrid);
% CB_pre = CB_factor * capacitance(xgrid);

capacitance_apical = struct( ...
    'dep', {{{'xpos'}}}, ...
    'dep_fun', {{{@(x) CA_factor * capacitance(x)}}} ...
);

capacitance_basal = struct( ...
    'dep', {{{'xpos'}}}, ...
    'dep_fun', {{{@(x) CB_factor * capacitance(x)}}} ...
);

%%

if args.plotflag
    
    xx = linspace(0,1,256);
    
    hfig = figure;
    hold on

    plot(xx, capacitance(xx), 'Color', 'b', 'LineStyle', '-')
    plot(xx, CA_factor * capacitance(xx), 'Color', 'b', 'LineStyle', '--')
    plot(xx, CB_factor * capacitance(xx), 'Color', 'b', 'LineStyle', ':')

    fig = Johnson_2011_data();
    fig(7).A.datasets(1).plotopt = [fig(7).A.datasets(1).plotopt, 'Color', 'r'];
    fig(7).A.datasets(2).plotopt = [fig(7).A.datasets(2).plotopt, 'Color', 'r'];
    leg = Johnson_2011_plot(fig(7), 'A', char_position);

    data = IHC.Johnson_2015.capacitance;
    errorbar(data.x, data.C, data.yerr, ...
        'Color', 'b', 'LineStyle', 'none', 'Marker', 's');

    set(gca,'XScale', 'linear')

    legend([{'IHC', 'IHC apical', 'IHC basolateral'}, leg, {'gerbil'}], ...
           'Location', 'NorthWest')
    title('IHC Capacitance');
    xlabel('relative distance from stapes');
    ylabel('C [pF]');

    if args.save_figures
        mySaveFig(hfig, 'IHC_capacitance', [], 'Results/HairCells/', default_tikz_options{:});
    end
else
    hfig = [];
end


%%

    function [b, a] = linear_regression(X, y)
        
        beta = X\y;
        
        a = beta(2);
        b = beta(1);
    end



end

