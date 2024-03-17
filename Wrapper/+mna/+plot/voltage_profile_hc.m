function [ ] = voltage_profile_hc( xx, fVolt, mnaopt, plotopt, runopt, analysisStartFrame)
%VOLTAGE_PROFILE_HC

sim_to_mV = Unit.conversionConstant(mnaopt.simulation_units.voltage, 'mV');

% ---------------------------------------------------------------------
% IHC
if plotopt.subplot == true
    hfig_ss_IHC_OHC = plotopt.figure;
    subplot(1,2,1);
else
    hfig_ss_IHC = plotopt.figure;
end

plot_steady_state_and_profile(xx, sim_to_mV * getVoltage(fVolt, 'IHC', mnaopt, 'transformed'), analysisStartFrame);

xlabel('BM position');
ylabel('Voltage [mV]');
title('IHC');

jj = find(strcmp({runopt.analysis_results.oc_mna.fullid}, 'IHCabsolute'), 1);
if isempty(jj)
    error('Expected to find a file with variable fullid = %s', 'IHCabsolute')
end
stats = matfile(runopt.analysis_results.oc_mna(jj).analysis_file, 'Writable', false);

VAR = sim_to_mV * getVoltage(fVolt, 'IHC', mnaopt, 'transformed');

hold on
if isfield('argmax_std', stats.oscillation_max)
    ind = stats.oscillation_max.argmax_std;
    plot(stats.x(:, ind), max(VAR(:,ind)), 'rx')
end
if isfield('char_frequency_model', stats.oscillation_max)
    char_frequency_model
    ind = stats.oscillation_max.char_frequency_model;
    plot(stats.x(:, ind), max(VAR(:,ind)), 'ro')
end
if isfield('max_delta_restricted', stats.oscillation_max)
    char_frequency_model
    ind = stats.oscillation_max.max_delta_restricted;
    plot(stats.x(:, ind), max(VAR(:,ind)), 'bx')
end
plot(stats.char_position, max(VAR(:,stats.char_position_ind)), 'bo')

% ---------------------------------------------------------------------
% OHC
if plotopt.subplot == true
    subplot(1,2,2);
else
    hfig_ss_OHC = plotopt.figure;
end

plot_steady_state_and_profile(xx, ...
    sim_to_mV * getVoltage(fVolt, 'OHC', mnaopt, 'transformed'),analysisStartFrame);

xlabel('BM position');
ylabel('Voltage [mV]');
title('OHC');

% ---------------------------------------------------------------------
% save figures
if runopt.save_figures        
        PLOTARG = {runopt.path.oc_mna, runopt.path.oc_mna, ...
            'height', '0.2\textwidth', ...
            'width', '0.4\textwidth', ...
            'showInfo', false};

    if plotopt.subplot == true
        mySaveFig(hfig_ss_IHC_OHC, 'ss_volt_IHC_OHC', PLOTARG{:});
    else
        mySaveFig(hfig_ss_IHC, 'ss_volt_IHC', PLOTARG{:});
        mySaveFig(hfig_ss_OHC, 'ss_volt_OHC', PLOTARG{:});
    end
end

end

