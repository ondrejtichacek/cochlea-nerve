function [ capacitance_apical, capacitance_basal, hfig ] = Capacitance( ...
    xgrid, simulation_units, args)
%CAPACITANCE
%OHC

arguments
    xgrid (:,1) double = linspace(0, 1, 300)'
    simulation_units (1,1) struct = struct()
    args.description = 'nonlinear'
    args.plotflag = false
    args.save_figures = false
end

%%

% [char_frequency, char_position] = frequency_map();
[char_position, char_frequency] = characteristic_frequency_model('Li');

mV_to_sim = Unit.conversionConstant('mV', simulation_units.voltage);
sim_to_mV = Unit.conversionConstant(simulation_units.voltage, 'mV');

default_tikz_options = { ...
    'relativeDataPath', 'img/HairCells/', ...
    'extraaxisoptions', {'legend style={font=\tiny}'}, ...
    'parseStrings', true};

%%

fig = Johnson_2011_data();

ref_V = -85e-3; % voltage, at which Johnson2011 measured the capacitance;

[sorted_cf, sorted_cap] = combineDatasets(fig(7).A.datasets, 'rat', 'gerbil');

FCi = griddedInterpolant(log10(sorted_cf), sorted_cap, 'pchip', 'linear');
P = polyfit(log10(sorted_cf),sorted_cap,1);
FCl = @(cf) P(1)*log10(cf) + P(2);

[~, ~, cf, cap] = combineDatasets(fig(7).A.datasets, 'rat', 'gerbil');

% extend the data manually to assure 'good' fit
extra_cf = [30, 25000];

extra_cap = [19, 4.5];

x0C = [18.6484  -14.7955    3.5674    0.2263];

[FCb, argss] = fitting.boltzmann([extra_cf, cf(:)'], [extra_cap, cap(:)'], x0C, optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4), @log10);

if args.plotflag
    c = logspace(log10(50), log10(20000), 100);
    
    hfig = figure;
    hold on

    legendtext = Johnson_2011_plot(fig(7), 'A');

    plot(extra_cf, extra_cap, 'Color', 'k', 'LineStyle', 'none', 'Marker', 'x')
    
    plot(c, FCi(log10(c)), 'Color', 'k', 'LineStyle', ':')
    plot(c, FCl(c), 'Color', 'k', 'LineStyle', '--')
    plot(c, FCb(c), 'Color', 'k', 'LineStyle', '-.')

    legend([legendtext, {'manual data', 'pchip interpolation', 'linear fit', 'boltzmann fit'}], 'Location', 'SouthWest') 
    
    title('OHC capacitance')
    
    if args.save_figures
        mySaveFig(hfig, 'OHC_capacitance', [], 'Results/HairCells/', default_tikz_options{:});
    end
end

OHC.capacitance_fit.interp = @(cf) FCi(log10(cf));
OHC.capacitance_fit.linear = FCl;
OHC.capacitance_fit.boltzmann = FCb;

%% Constant capacitance

% Johnson 2011 claims, that the ratio CA/CB of apical and basilar
% capacitance is r = 0.2 uniformly across the BM
% CA = CB*r, CB = CA/r

capacitance_ratio = 0.2;

% Also, we know that the total membrane capacitance CM = CA + CB
% Therefore
% CA = CM / (1 + 1/r);
% CB = CM / (1 + r);

CA_factor = 1/(1 + 1/capacitance_ratio);
CB_factor = 1/(1 + capacitance_ratio);

% The OHC capacitance is not uniform along the tonotopy axis

% extract the function to prevent hashing error
capacitance_tmp = OHC.capacitance_fit.boltzmann;
% capacitance_tmp = OHC.capacitance_fit.linear;
% capacitance_tmp = OHC.capacitance_fit.interp;

capacitance = @(x) capacitance_tmp(char_frequency(x));

CA_pre = CA_factor * capacitance(xgrid);
CB_pre = CB_factor * capacitance(xgrid);

%% Extenson -- Nonlinear capacitance

% DOI: 10.1002/cm.20423

Q_max = 978.9; % +- 182.6 fC
alpha = 1/31.1; % +- 4.9 mV^-1
Vhalf = -65.2; % +- 29.4 mV
Clin = 5.69; % +- 0.83 pF
Cnl = 7.9; % +- 1.6 pF
z = 0.79; % +- 0.12

Cnl_times_4 = Q_max * alpha; % now in pF

alpha = alpha / mV_to_sim;
Vhalf = Vhalf * mV_to_sim;


NLC_max = NLC(Vhalf, alpha, Vhalf, Cnl_times_4) + Clin; 

assert(abs(NLC_max - Cnl_times_4 / 4 - Clin)/Clin < 1e-6)

NLC_rel = @(V) (NLC(V, alpha, Vhalf, Cnl_times_4) + Clin) / NLC_max;


if args.plotflag
    V = -150:100;
    
    hfig = figure;
    plot(V, 100*NLC_rel(V* mV_to_sim), 'k');
    xlabel('Membrane potential (mV)');
    ylabel('Relative membrane capacitance (%)');
else
    hfig = [];
end

%%
% https://doi.org/10.1038/s41598-020-63201-6
% Fig 2
% C_m (peak) @ 18.89 kHz = 52.60 fF (patch)
% C_m (peak) @ 0.659 kHz = 133.60 fF (patch)

xxx = [char_position(18890); char_position(659)];
yyy = [52.6/133.6; 1];
X = [ones(length(xxx),1), xxx];
b = X\yyy;
NLC_freq = @(x) b(1) + b(2)*x;

assert(abs(NLC_max - Cnl_times_4 / 4 - Clin)/Clin < 1e-6)

NLC_rel_v2 = @(x,V) (NLC(V, alpha, Vhalf, Cnl_times_4) .* NLC_freq(x)  + Clin) / NLC_max;

NLC_freq_x_pre = NLC_freq(xgrid);

NLC_rel_v2_pre = @(V) (NLC(V, alpha, Vhalf, Cnl_times_4) .* NLC_freq_x_pre  + Clin) / NLC_max;


if args.plotflag
    V = -150:100;
    
    hfig = figure;
    hold on
    xx = linspace(0,1,10);
    cmp = flip(rwb(numel(xx)));
    for i = 1:numel(xx)
        x = xx(i);
        plot(V, 100*NLC_rel_v2(x, V* mV_to_sim), ...
            'Color', cmp(i,:), ...
            'LineStyle', '-');
    end

    xlabel('Membrane potential (mV)');
    ylabel('Relative membrane capacitance (%)');
    
    cbarTL = cellfun(@(x) sprintf('%d', round(Frequency(char_frequency(x), 'Hz').Hz)), num2cell(xx), 'UniformOutput', false);
            
    manual_colorbar(cmp, cbarTL, ...
        'orientation', 'vertical', ...
        'reversed', true, ...
        'label', sprintf('f (Hz)'));

else
    hfig = [];
end

%% 
switch args.description
    case 'linear'
        
        CA = @(x) CA_factor * capacitance(x);
        CB = @(x) CB_factor * capacitance(x);
        
        capacitance_apical = struct( ...
            'dep', {{{'xpos'}}}, ...
            'dep_fun', {{{CA}}} ...
        );
        
        capacitance_basal = struct( ...
            'dep', {{{'xpos'}}}, ...
            'dep_fun', {{{CB}}} ...
        );
        
        
    case 'nonlinear'
        
        fac = 1 / NLC_rel(ref_V);

        CA = @(x) fac * CA_factor * capacitance(x);
        CB = @(x,V) fac * NLC_rel(V) .* CB_factor .* capacitance(x);
        CBb = @(V) fac * NLC_rel(V) .* CB_pre;
        
        capacitance_apical = struct( ...
            'dep', {{{'xpos'}}}, ...
            'dep_fun', {{{CA}}} ...
        );
        
        capacitance_basal = struct( ...
            'dep', {{{'xpos', 'vohc'}}}, ...
            ...
            'dep_fun', {{{CB}}}, ...
            ...
            'dep_fun_args', {{{{'xpos', 'vohc'}}}}, ...
            'pre', struct( ...
            'dep_fun', struct( ...
                ...
                'fun', CBb, ...
                ...
                'args', {{'vohc'}})) ...
        );
    
    case 'nonlinear_freq'
        
        fac = 1 / NLC_rel(ref_V);

        CA = @(x) fac * CA_factor * capacitance(x);
        CB = @(x,V) fac * NLC_rel_v2(x,V) .* CB_factor .* capacitance(x);
        CBb = @(V) fac * NLC_rel_v2_pre(V) .* CB_pre;
        
        capacitance_apical = struct( ...
            'dep', {{{'xpos'}}}, ...
            'dep_fun', {{{CA}}} ...
        );
        
        capacitance_basal = struct( ...
            'dep', {{{'xpos', 'vohc'}}}, ...
            ...
            'dep_fun', {{{CB}}}, ...
            ...
            'dep_fun_args', {{{{'xpos', 'vohc'}}}}, ...
            'pre', struct( ...
            'dep_fun', struct( ...
                ...
                'fun', CBb, ...
                ...
                'args', {{'vohc'}})) ...
        );

    case 'nonlinear_freq_ss'
        
        fac = 1 / NLC_rel(ref_V);

        CA = @(x) fac * CA_factor * capacitance(x);
        CB = @(x,V) fac * NLC_rel_v2(x,V) .* CB_factor .* capacitance(x);
        CBb = @(V) fac * NLC_rel_v2_pre(V) .* CB_pre;
        
        capacitance_apical = struct( ...
            'dep', {{{'xpos'}}}, ...
            'dep_fun', {{{CA}}} ...
        );
        
        capacitance_basal = struct( ...
            'dep', {{{'xpos', 'vohc_ss'}}}, ...
            ...
            'dep_fun', {{{CB}}}, ...
            ...
            'dep_fun_args', {{{{'xpos', 'vohc_ss'}}}}, ...
            'pre', struct( ...
            'dep_fun', struct( ...
                ...
                'fun', CBb, ...
                ...
                'args', {{'vohc_ss'}})) ...
        );
    
    case 'nonlinear_freq_ss_v2'
        
        fac = 1 / NLC_rel(ref_V);

        CA = @(x) fac * CA_factor * capacitance(x);        
        CB = @(x,V) fac * NLC_rel_v2(x,V) .* CB_factor .* capacitance(x);
        CBb = @(V) fac * NLC_rel_v2_pre(V) .* CB_pre;
        
        ss_voltage = load('ss_voltage/yxj6cwdmzoqu5pkoqgimwv73sxnbmkila.mat');

        CB_eval = CBb(ss_voltage.vohc_ss);
        CB_fun = griddedInterpolant(xgrid, CB_eval);

        capacitance_apical = struct( ...
            'dep', {{{'xpos'}}}, ...
            'dep_fun', {{{CA}}} ...
        );        

        capacitance_basal = struct( ...
            'dep', {{{'xpos'}}}, ...
            'dep_fun', {{{@(x) CB_fun(x)}}}, ...
            'pre', struct( ...
                    'dep_fun', struct( ...
                        'fun', CB_eval, ...
                        'args', {{}})) ...
        );

    otherwise
        error(args.description)
end

end
function cap = NLC(V, alpha, Vhalf, Cnl_times_4)
    % NLC = @(V) Cnl_times_4 ./ (exp(alpha * (V - Vhalf)) .* ((1 + exp(-alpha * (V - Vhalf))).^2)) + Clin;
    Z = alpha * (V - Vhalf);
    cap = Cnl_times_4 ./ (exp(Z) .* ((1 + exp(-Z)).^2));
end

