function [ ] = current_profile( xx, fCurr, mnaopt, plotopt, runopt, analysisStartFrame, relative )
%CURRENT_PROFILE 

if plotopt.subplot == true
    f1 = plotopt.figure;
end

% ---------------------------------------------------------------------
% IHC
if plotopt.subplot == true
    subplot(2,2,1);
else
    sf1 = plotopt.figure;
end

SCALING = 1e12; % from A to pA
SCALING = SCALING / mnaopt.NumstacksScaling; % approximate single IHC 

plot_steady_state_and_profile(xx, fCurr.IHC, analysisStartFrame, SCALING, relative);

xlabel('BM position');
ylabel('Current [pA]');
title('single IHC');

% IHC_MET
if plotopt.subplot == true
    subplot(2,2,2);
else
    sf1 = plotopt.figure;
end

SCALING = 1e12; % from A to pA
SCALING = SCALING / mnaopt.NumstacksScaling; % approximate single IHC 

plot_steady_state_and_profile(xx, fCurr.IHC_MET, analysisStartFrame, SCALING, relative);

xlabel('BM position');
ylabel('Current [pA]');
title('single IHC MET');

% ---------------------------------------------------------------------
% OHC
if plotopt.subplot == true
    subplot(2,2,3)
else
    sf2 = plotopt.figure;
end

SCALING = 1e12; % from A to pA
SCALING = SCALING / mnaopt.NumstacksScaling; % approximate single OHC 
SCALING = SCALING / mnaopt.num_ohc_per_cs; % approximate single OHC 

plot_steady_state_and_profile(xx, fCurr.OHC, analysisStartFrame, SCALING, relative);

xlabel('BM position');
ylabel('Current [pA]');
title('single OHC');

% OHC_MET
if plotopt.subplot == true
    subplot(2,2,4)
else
    sf2 = plotopt.figure;
end

SCALING = 1e12; % from A to pA
SCALING = SCALING / mnaopt.NumstacksScaling; % approximate single OHC 
SCALING = SCALING / mnaopt.num_ohc_per_cs; % approximate single OHC 

plot_steady_state_and_profile(xx, fCurr.OHC_MET, analysisStartFrame, SCALING, relative);

xlabel('BM position');
ylabel('Current [pA]');
title('single OHC MET');

% ---------------------------------------------------------------------
% save figures
if runopt.save_figures

    PLOTARG = {runopt.path.oc_mna, runopt.path.oc_mna, ...
            'height', '0.2\textwidth', ...
            'width', '0.4\textwidth', ...
            'showInfo', false};

    if plotopt.subplot == true
        mySaveFig(f1, 'ss_current_IHC_OHC', PLOTARG{:});
    else

        mySaveFig(sf1, 'ss_current_IHC', PLOTARG{:});
        mySaveFig(sf2, 'ss_current_OHC', PLOTARG{:});

    end
end

end

