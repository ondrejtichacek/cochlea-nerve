function [ ] = voltage_profile( xx, fVolt, mnaopt, plotopt, runopt, analysisStartFrame )
%VOLTAGE_PROFILE

sim_to_mV = Unit.conversionConstant(mnaopt.simulation_units.voltage, 'mV');


if plotopt.subplot == true
    f2 = plotopt.figure;
end

% ---------------------------------------------------------------------
% Scala Media
if plotopt.subplot == true
    subplot(2,3,3);
else
    sf3 = plotopt.figure;
end

plot_steady_state_and_profile(xx, sim_to_mV * getVoltage(fVolt, 'ScalaMedia', mnaopt, 'transformed'), analysisStartFrame);

xlabel('BM position');
ylabel('Voltage [mV]');
title('Scala Media');

% ---------------------------------------------------------------------
% IHC intracellular
if plotopt.subplot == true
    subplot(2,3,1);
else
    sf4 = plotopt.figure;
end

plot_steady_state_and_profile(xx, sim_to_mV * getVoltage(fVolt, 'IHC_ic', mnaopt, 'transformed'), analysisStartFrame);

xlabel('BM position');    
ylabel('Voltage [mV]');
title('IHC intracellular');

% ---------------------------------------------------------------------
% OHC intracellular
if plotopt.subplot == true
    subplot(2,3,2);
else
    sf5 = plotopt.figure;
end

plot_steady_state_and_profile(xx, sim_to_mV * getVoltage(fVolt, 'OHC_ic', mnaopt, 'transformed'), analysisStartFrame);

xlabel('BM position');
ylabel('Voltage [mV]');
title('OHC intracellular');

% ---------------------------------------------------------------------
% IHC extracellular
if plotopt.subplot == true
    subplot(2,3,4);
else
    sf6 = plotopt.figure;
end

plot_steady_state_and_profile(xx, sim_to_mV * getVoltage(fVolt, 'IHC_ec', mnaopt, 'transformed'), analysisStartFrame);

xlabel('BM position');
ylabel('Voltage [mV]');
title('IHC extracellular');

% ---------------------------------------------------------------------
% OHC extracellular
if plotopt.subplot == true
    subplot(2,3,5);
else
    sf7 = plotopt.figure;
end

plot_steady_state_and_profile(xx, sim_to_mV * getVoltage(fVolt, 'OHC_ec', mnaopt, 'transformed'), analysisStartFrame);

xlabel('BM position');
ylabel('Voltage [mV]');
title('OHC extracellular');

% ---------------------------------------------------------------------
% Spiral Ligament
if plotopt.subplot == true
    subplot(2,3,6);
else
    sf8 = plotopt.figure;
end

plot_steady_state_and_profile(xx, sim_to_mV * getVoltage(fVolt, 'SpiralLigament', mnaopt, 'transformed'), analysisStartFrame);

xlabel('BM position');
ylabel('Voltage [mV]');
title('Spiral Ligament');

% ---------------------------------------------------------------------
% save figures
if runopt.save_figures

    PLOTARG = {runopt.path.oc_mna, runopt.path.oc_mna, ...
        'height', '0.2\textwidth', ...
        'width', '0.4\textwidth', ...
        'showInfo', false};

    if plotopt.subplot == true
        mySaveFig(f2, 'ss_volt', PLOTARG{:});
    else
        mySaveFig(sf3, 'ss_volt_SM', PLOTARG{:});
        mySaveFig(sf4, 'ss_volt_IHC_ic', PLOTARG{:});
        mySaveFig(sf5, 'ss_volt_OHC_ic', PLOTARG{:});
        mySaveFig(sf6, 'ss_volt_IHC_ec', PLOTARG{:});
        mySaveFig(sf7, 'ss_volt_OHC_ec', PLOTARG{:});
        mySaveFig(sf8, 'ss_volt_SL', PLOTARG{:});
    end
end


end

