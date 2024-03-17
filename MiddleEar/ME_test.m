function [analysis_variables, CACHE] = ME_test(varargin)
arguments (Repeating)
    varargin
end
arguments           
    ...args.tf_extra (1,1) Time ... % extra time after defult analysis start time
    ...    = Time(40, 'ms')
end

args = struct();

[args, opt, memopt, paropt] = common_opts(args, varargin{:});

%% Plotting options

plotopt = plotOpt('MNA', ...
            ...
            'default_plt_action', false, ...
            ...
            ...'subplot', 0);
            'subplot', true);
        
% plotopt.hide_figures_while_plotting();

%% Default options

[~, stimulus] = devopts.defaults(args, 'multiple', plotopt);

%% Select units for simulation

% units
% si: F * V / s + S * V = A
%     uF * mV / ms + mS * mV = uA

% si: H * A / s + V = V
%     H * uA / ms + mV = mV

simulation_units = struct( ...
    'time', 's', ...
    'voltage', 'V', ...
    'current', 'A', ...
    'conductance', 'S', ...
    'inductance', 'H', ...
    'capacitance', 'F');

% simulation_units = struct( ...
%     'time', 'ms', ...
%     'voltage', 'mV', ...
%     'current', 'uA', ...
%     'conductance', 'mS', ...
%     'inductance', 'H', ...
%     'capacitance', 'uF');

%% Middle Ear options

updated_parameters = struct( ...
    'Rcm', 0.002937171395014418, ...
    'La', 0.007821011726873243, ...
    'Ra',   0.1781859541682067, ...
    'Rti1', 0.005717623809642999, ...
    'Cti1', 3.665604886531633e-06, ...
    'Cti2', 2.331753061471045e-07, ...
    'Rti2',  0.01816201912960247, ...
    'Rti3', 0.0001154419087022123, ...
    'Lti',    0.015563307309665, ...
    'Cti3', 5.180209724606197e-07, ...
    'Cte', 9.019459088027458e-07, ...
    'Lte',  0.04823092240082723, ...
    'Rte',  0.02716531759951024, ...
    'Cis', 4.61959539247965e-08, ...
    'Ris',  0.04427903532306574, ...
    'Cst', 1.823021558987398e-06, ...
    'Rla',   0.0162895857659625, ...
    'Cla', 3.991813402702068e-07, ...
    'Rh', 0.001136017714952076, ...
    'Lh',  0.02150377270840547, ...
    'Ccp', 4.694955211672143e-06 ...
);

circuit = ME_circuit_PBLL(simulation_units, ...
    'UpdatedParameters', updated_parameters); IND = 17;

% circuit = ME_circuit_PBLL(simulation_units); IND = 17;
% circuit = ME_circuit_Zwislocki(simulation_units); IND = 11;

midopt = midOpt( ...
    'pv_TM', 1e6, ...
    'circuit', circuit, ...
    'simulation_units', simulation_units, ...
    'IND', IND);


%% Run options
runopt = runOpt( ...
    'save_figures', true, ...
    'figureSaveDir', fullfile(opt.cochleadir, 'Results', 'ME'), ...
    'CodeVersion', codeVersion(), ...
    'Debug', false, ...
    'draw_gui', false, ...
    'analysisStart', Time(0), ...
    'waitbarFunctionAvailable', false, ... % disable waitbar
    'verbose', 2 );
% verbose levels:   0 ... just matlab errors & warnings
%                   1 ... low verbosity
%                   2 ... normal verbosity
%                   3 ... high verbosity


% default action for do is TRUE
skip_do = { ...
    ...'mid', ...
    };

% default action for recalculate is FALSE
do_recalculate = { ...
    ...'mid', ...
    };

% default action for plot is FALSE
do_plot = { ...
    ...'mid', ...
    };

% update in runopt
runopt.update_struct('do', skip_do, false);
runopt.update_struct('recalculate', do_recalculate, true);
runopt.update_struct('plot', do_plot, true);


runopt.save_figures = false;
% runopt.save_figures = true;

export_flag = false;
% export_flag = true;

%% Plotting options
plotopt = plotOpt('ME');

%%

% AMPLITUDE = 0:5:120;
AMPLITUDE = 60;

FREQUENCY = round(erbspace(100,20000,96));
FREQUENCY = FREQUENCY(1:4:end);

[Configurations, Conf_nums] = ParameterProduct( ...
    'Amplitude', AMPLITUDE, 'Frequency', FREQUENCY);
SIM_OPTS = {midopt};

stats = batchGathering(@MiddleEar_MNA_statistics, 1, Configurations, ...
            stimulus, SIM_OPTS, ...
            runopt, opt, memopt, paropt);

       
res = @(x) reshape(x, numel(AMPLITUDE), numel(FREQUENCY));
       
% Y = zeros(numel(AMPLITUDE), numel(FREQUENCY));
% for i = 1:length(AMPLITUDE)
%     for j = 1:length(FREQUENCY)        
%         ind = sub2ind(Conf_nums,i,j);
%         Y(i,j) = stats{ind}.max;
%     end
% end

stats = [stats{:}];

Y = res([stats.max]);

phase_lag_fft = res([stats.phase_lag_fft]);
phase_lag_hilbert = res([stats.phase_lag_hilbert]);

phase_lag = phase_lag_fft;
% phase_lag = phase_lag_hilbert;

%% Transfer function

exp_data = {'Aibara_2001', 'Chien_2009'};
% exp_data = {'Chien_2009'};

hfig = plotStapesVelocityTransferFunction( ...
    'do_plot_selected', exp_data);

figure(hfig(1));

for i = 1:numel(AMPLITUDE)
    plot(FREQUENCY,mag2db(Y(i,:)/max(Y(i,:))), ...
        'Color', 'k', 'LineWidth', 1, 'DisplayName','model');
end

ylabel('Velocity (dB)')

% Save figure

FIGNAME = sprintf('ME-filtering-PBLL-dev');

PATH = 'img/ME/';

if runopt.save_figures
    analysis.plot.export_figure(hfig(1), FIGNAME, PATH, export_flag, runopt)
end


figure(hfig(2));

for i = 1:numel(AMPLITUDE)
%     plot(FREQUENCY, unwrap(phase_lag_fft(i,:)));
    plot(FREQUENCY,-180 + 180/pi*unwrap(phase_lag(i,:)), ...
        'Color', 'k', 'LineWidth', 1, 'DisplayName','model');
%     plot(FREQUENCY,180 + 180/pi*unwrap(phase_lag(i,:)), ...
%         'Color', 'k', 'LineWidth', 2, 'DisplayName','model');
end

legend('off')

% Save figure

FIGNAME = sprintf('ME-phase-lag-PBLL-dev');

PATH = 'img/ME/';

if runopt.save_figures
    analysis.plot.export_figure(hfig(2), FIGNAME, PATH, export_flag, runopt)
end


%%

% AMPLITUDE = 0:5:120;
AMPLITUDE = 0:40:120;

% FREQUENCY = round(erbspace(100,20000,48));
FREQUENCY = round(erbspace(100,20000,48));
FREQUENCY = FREQUENCY(1:9:end);
% FREQUENCY = 2000;

[Configurations, Conf_nums] = ParameterProduct( ...
    'Amplitude', AMPLITUDE, 'Frequency', FREQUENCY);
SIM_OPTS = {midopt};

stats = batchGathering(@MiddleEar_MNA_statistics, 1, Configurations, ...
            stimulus, SIM_OPTS, ...
            runopt, opt, memopt, paropt);

res = @(x) reshape(x, numel(AMPLITUDE), numel(FREQUENCY));
       
% Y = zeros(numel(AMPLITUDE), numel(FREQUENCY));
% for i = 1:length(AMPLITUDE)
%     for j = 1:length(FREQUENCY)        
%         ind = sub2ind(Conf_nums,i,j);
%         Y(i,j) = stats{ind}.max;
%     end
% end

stats = [stats{:}];

Y = res([stats.max]);

%%

hfig = figure;
hold on
cmap = rwb(numel(FREQUENCY));
for i = 1:numel(FREQUENCY)
    plot(AMPLITUDE,Y(:,i), 'Color', cmap(i,:));
end
set(gca, 'YScale', 'log')
xlabel('Amplitude (dB SPL)');
ylabel('$\max(V_{out})$ [V]');
legendtext = arrayfun(@(x) sprintf('%d Hz', x), FREQUENCY, 'UniformOutput', false);
legend(legendtext{:}, 'Location', 'SouthEast')
    

% Save figure

FIGNAME = sprintf('ME-magnitude-PBLL');

PATH = 'img/ME/';

if runopt.save_figures
    analysis.plot.export_figure(hfig, FIGNAME, PATH, export_flag, runopt)
end

end