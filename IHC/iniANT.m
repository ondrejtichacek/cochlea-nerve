function [ potential, displacement, TimeVoltage, Signal, TimeSignal, antopt, runopt ] = iniANT( ...
        stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt )
%INIANT Initialization of ANT


% warning('off', 'MATLAB:hg:NoDisplayNoFigureSupportSeeReleaseNotes') % supress waitbar warning in nodisplay mode

% Reconstruct results
% [resFiles, runoptMNA] = MNA( stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt );

results = MNA( ...
    stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, ...
    'all_analysis_files_must_exist', ...
    false ...
    ...true ...
    );

BMx = results.result_files.fMech.BMx;

VIHC = Voltage(getVoltage(results.result_files.fVolt, 'IHC', mnaopt, 'transformed'), mnaopt.simulation_units.voltage);

VIHC = VIHC.mV;

TimeVoltage = Time(results.result_files.fT.T, results.result_files.fT.unit);

if false % old version
    % fn = 'fBM';
    fn = 'fMech';

    Signal = results.result_files.(fn).signal;
    TimeSignal = Time(results.result_files.(fn).time, results.result_files.(fn).time_unit);
else
    TimeSignal = TimeVoltage;
    Signal = stimulus.eval(TimeSignal);
end
    
% Update loaded 'runopt' (from MNA) with the input 'runopt'
% runoptMNA.updateWithObject(runopt);
% runopt = runoptMNA;

disp_title( 'Stored cochlear simulation parameters' );
% displayStats( runopt.verbose, fullfile( runopt.path.oc_mna, 'options.mat') )


% Cross-sections to excite
sig_test = BMx;
% sig_test = VIHC;

[s, x] = slicesToExcite( antopt.slices, mnaopt.Numstacks, sig_test, results, mnaopt.xgrid );
antopt.update('slicesToExcite', s);
antopt.update('positionsToExcite', x);

potential = VIHC(:,antopt.slicesToExcite);
displacement = BMx(:,antopt.slicesToExcite);

end

