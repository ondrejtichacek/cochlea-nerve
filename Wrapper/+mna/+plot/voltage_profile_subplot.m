function [ ] = voltage_profile_subplot( xx, volt_fun, variables, mnaopt, plotopt, runopt, analysisStartFrame )
%VOLTAGE_PROFILE_SUBPLOT 

if plotopt.subplot == true
    f2 = plotopt.figure;
end

for i = 1:numel(variables)
    if plotopt.subplot == true
        subplot(2,ceil(numel(variables)/2),i);
    else
        sf3 = plotopt.figure;
    end

    plot_steady_state_and_profile(xx, volt_fun(i), analysisStartFrame);

    xlabel('BM position');
    ylabel('Voltage [mV]');
    title(variables{i});

end

end

