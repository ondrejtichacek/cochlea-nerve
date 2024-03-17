function [status] = MNAplotting( results, plotopt, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt )
%PLOTTING
arguments
    results (1,1) struct
    plotopt (1,1) plotOpt
    stimulus
    topt (1,1)
    midopt (1,1) midOpt
    mechopt (1,1) mechOpt
    mnaopt (1,1) mnaOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
end

doplot = plotopt.do.oc_electric;

%% Stimulus
if doplot.stimulus

    hfig = plotopt.figure;

    % try
        plot(stimulus);
    % catch ME
    %     warning('Stimulus plotting failed:')
    %     disp(ME)
    % end

    if runopt.save_figures
        plotopt.save_figure(hfig,  runopt.figureSaveDir, 'stimulus');
    end
end

%%

status = true;

resFiles = results.oc_mna;

% matfile alias
fT    = resFiles.fT;
fVolt = resFiles.fVolt;
fChannels = resFiles.fChannels;
fCurr = resFiles.fCurr;
if strcmp(mechopt.integration, 'electric')
    fMech = resFiles.fMech;
end

%% Plotting results

disp_title( sprintf('plotting MNA run') );
if runopt.save_figures
    if isempty(runopt.figureSaveDir)
        runopt.figureSaveDir = runopt.path.oc_mna;
    end
    fprintf('Figures will be saved to:\n%s\n', runopt.figureSaveDir);
end

%%

sim_to_mV = Unit.conversionConstant(mnaopt.simulation_units.voltage, 'mV');
mV_to_sim = Unit.conversionConstant('mV', mnaopt.simulation_units.voltage);

% starting time-index of the analysis
t1 = find(runopt.analysisStart <= Time(fT.T, fT.unit), 1);

%% plotting options

switch plotopt.xAxisValues
    case 'stacks'
        xx = 1:mnaopt.Numstacks;
    case 'x'
%         xx = fmech.x
%         xx = fBM.x;
        xx = mnaopt.xgrid;
end

% tmp fix
xx = xx(:);

SamplingFrequency = 'optimal';

%% 
if doplot.gui
    mech.gui_v2.loader(results, stimulus, mechopt, mnaopt);
end

%% clear stored figure handles

plotopt.figure_handles = [];

%% MNA settings

% =========================================================================
% CIRCUIT
% -------------------------------------------------------------------------
if doplot.MNAsettings    
    VoltIHC_0 = subsampleVoltage(fVolt, 'IHC', 1, mnaopt, plotopt, 'voltage');
    VoltOHC_0 = subsampleVoltage(fVolt, 'OHC', 1, mnaopt, plotopt, 'voltage');
    
    VoltIHC_0 = VoltIHC_0 * mV_to_sim;
    VoltOHC_0 = VoltOHC_0 * mV_to_sim;
    
    mna.plot.resistors(xx, VoltIHC_0, VoltOHC_0, mnaopt, plotopt, runopt, 'siemens');
    mna.plot.resistors(xx, VoltIHC_0, VoltOHC_0, mnaopt, plotopt, runopt, 'ohm');
	mna.plot.capacitors(xx, VoltIHC_0, VoltOHC_0, mnaopt, plotopt, runopt);
end
% =========================================================================

%% Maximal cross-section
V = sim_to_mV * getVoltage(fVolt, 'IHC', mnaopt, 'transformed'); % all voltages
V = V(t1:end,:);                                                 % data after analysis start
V_max_allpos = max(V,[],1);                                      % vector of maximum values at each BM position.
[~,ii_max] = max(V_max_allpos);                                  % index of the BM segment with Vmax
IHC_V = V(:,ii_max);                                             % temporal evolution of voltages at the BM segment with Vmax

%% Fourier transform

if doplot.FourierTransform
    load(fullfile(runopt.path.mech, 'BM.mat'), 'signal');

    if all([stimulus.samplingFrequency] == stimulus(1).samplingFrequency)
        fs = stimulus(1).samplingFrequency;
    else
        error('Stimulus sampling frequency inconsistent! %g\n',stimulus.samplingFrequency)
    end

    NN = length(signal); % Length of signal
    XX = ( 1:NN ) / fs;

    %     SampleTime = 1/44100;                     % Sample time
    %     t = (0:L-1)*SampleTime;                % Time vector
    [f,YY] = plotSignalFT( signal, fs );
end
[fV,VV] = plotSignalFT( IHC_V - mean(IHC_V), mnaopt.samplingFrequency );

if doplot.FourierTransform
    
    if plotopt.subplot == true
        f4 = plotopt.figure;
    else
        sf10 = plotopt.figure;
    end
    plot(XX.ms, signal(1:NN));
    title('Input sinal')
    xlabel('Time (sec)')
    xlabel('Time, [ms]');
    xlim([0 XX(end).ms]);
    
    
    if plotopt.subplot == true
        f5 = plotopt.figure;
        subplot(2,1,1);
    else
        sf11 = plotopt.figure;
    end
    
    
    % Plot single-sided amplitude spectrum.
    plotyy(f.Hz,YY,fV.Hz,VV)
    xlabel('Frequency (Hz)')
    
    legend({'FT of input signal','FT of IHC voltage signal'}, 'location', 'best');
    
    if plotopt.subplot == true
        subplot(2,1,2);
    else
        sf12 = plotopt.figure;
    end
    
    [~,xxx] = max(YY);
    d = 0.04*f(xxx).Hz+eps;
    
    ax = plotyy(f.Hz,YY,fV.Hz,VV);
    set(ax, 'xlim', [f(xxx).Hz - d , f(xxx).Hz + d]);
    
    xlabel('Frequency (Hz)')
    
    %% save figures
    if runopt.save_figures
        if plotopt.subplot == true
            plotopt.save_figure(f4,  runopt.figureSaveDir, 'figure_4');
            plotopt.save_figure(f5,  runopt.figureSaveDir, 'figure_5');
        else
            plotopt.save_figure(sf10,  runopt.figureSaveDir, 'sfigure_10');
            plotopt.save_figure(sf11,  runopt.figureSaveDir, 'sfigure_11');
            plotopt.save_figure(sf12,  runopt.figureSaveDir, 'sfigure_12');
        end
    end
end

%% undersampling

[~,m] = max(VV);

fsampleminrequired = Frequency(2*fV(m).Hz, 'Hz');

fsampling = mnaopt.samplingFrequency;

% SamplingFrequency = 'minimal';
% SamplingFrequency = fsampling;
% SamplingFrequency = 44.1e3;
% SamplingFrequency = 'maximal';


switch SamplingFrequency
    case 'optimal'
        fsampldesired = exp(1)*fsampleminrequired/2; % desired sampling frequency is somewhat larger than the minimal sampling frequency according to the Shannon-Nyquist theorem
    case 'minimal'
        fsampldesired = fsampleminrequired;
    case 'maximal'
        fsampldesired = fsampling;
    otherwise
        fsampldesired = SamplingFrequency;
        
end

fsampldesired = Frequency(20, 'kHz');

if isa(stimulus, 'PureTone')
    fsampldesired = 16 * stimulus.frequency;
end

% if isempty(SamplingSkip)
%     SamplingSkip = max([1, floor(fsampling/fsampldesired)]);
% end
% 
% 
% originalSampligTime = TIME(2,1) - TIME(1,1);
% 
% ts = 1:SamplingSkip:size(TIME,1); % time samples
% 
% % ts = ts+50;
% % 
% % ts(TIME(ts) > 0.06 | TIME(ts) < 0.01) = [];
% 
% T = TIME(ts,:);
% 
% visualisedSampligTime = T(2,1) - T(1,1);
% 
% warning(...
%     ['ploted signal undersampled for visualisation by factor %g\n', ...
%     '\t original sampling time/frequency %g s / %g Hz\n', ...
%     '\t visualised sampling time/frequency %g s / %g Hz\n', ...
%     '\t minimal required sampling time/frequency %g s / %g Hz (Nyquist)\n'], ...
%     SamplingSkip, ...
%     originalSampligTime, 1/originalSampligTime, ...
%     visualisedSampligTime, 1/visualisedSampligTime, ...
%     1/fsampleminrequired, fsampleminrequired);

%% TIME

time = Time(fT.T, fT.unit);

[T, ts, t0, tf] = plotopt.time_samples(time, topt, fsampldesired);

ts = ts(ts>=t1);
t0 = time(t1);
T = T(T>=time(t1).ms);
time = time(t1:end,:);

[Xgrid,Tgrid] = meshgrid(xx, T);

[Xgrid_mech,Tgrid_mech] = meshgrid(xx, T);
% [Xgrid_mech,Tgrid_mech] = meshgrid(1:mechopt.Numstacks, T);

%% Maximal cross-section plot

if doplot.MaximalCrossection
    mna.plot.max_crossection(time, T, t0, tf, IHC_V, mnaopt, plotopt, runopt)
end


%% Current
if doplot.ApicalCurrentMesh
    CurrIHC = plotopt.subsample(fCurr, 'IHC_MET', ts, ':');
    CurrIHC = CurrIHC * Current(1, mnaopt.simulation_units.current).pA;
    CurrIHC = CurrIHC / mnaopt.NumstacksScaling; % approximate single IHC
    
    CurrOHC = plotopt.subsample(fCurr, 'OHC_MET', ts, ':');
    CurrOHC = CurrOHC * Current(1, mnaopt.simulation_units.current).pA;
    CurrOHC = CurrOHC / mnaopt.NumstacksScaling; % approximate single OHC
    CurrOHC = CurrOHC / mnaopt.num_ohc_per_cs; % approximate single OHC
    
    mna.plot.current_mesh_hc(Xgrid, Tgrid, CurrIHC, CurrOHC, t0, tf, plotopt, runopt, 'pA');
end
if doplot.BasalCurrentMesh
    CurrIHC = plotopt.subsample(fCurr, 'IHC', ts, ':');
    CurrIHC = CurrIHC * Current(1, mnaopt.simulation_units.current).pA;
    CurrIHC = CurrIHC / mnaopt.NumstacksScaling; % approximate single IHC
    
    CurrOHC = plotopt.subsample(fCurr, 'OHC', ts, ':');
    CurrOHC = CurrOHC * Current(1, mnaopt.simulation_units.current).pA;
    CurrOHC = CurrOHC / mnaopt.NumstacksScaling; % approximate single OHC
    CurrOHC = CurrOHC / mnaopt.num_ohc_per_cs; % approximate single OHC
    
    mna.plot.current_mesh_hc(Xgrid, Tgrid, CurrIHC, CurrOHC, t0, tf, plotopt, runopt, 'pA');
end

%% Current Steady state

if doplot.CurrentSteadyState
    mna.plot.current_profile(xx, fCurr, mnaopt, plotopt, runopt, t1, false);
    mna.plot.current_profile(xx, fCurr, mnaopt, plotopt, runopt, t1, true);
end


%% Voltage HC

    function plot_voltage_mesh(SEL, ss)
        V_ss = subsampleVoltage(fVolt, SEL, 1, mnaopt, plotopt, 'voltage');
        V_t = subsampleVoltage(fVolt, SEL, ts, mnaopt, plotopt, 'voltage');

        if ss == true
            V_t = V_t - V_ss;
            V0 = 0;
        else
            V0 = [];
        end

        hfig = mna.plot.voltage_mesh_one(Xgrid, Tgrid, V_t, t0, tf, plotopt, V0);
        
        title(SEL)
        
        % save figure        
        if runopt.save_figures
            fname = sprintf('Volt_%s', SEL);
            plotopt.save_figure(hfig, runopt.figureSaveDir, fname);
        end
    end

if doplot.VoltageMesh_IHC_rel
    plot_voltage_mesh('IHC', true)
    plot_voltage_mesh('IHC_MET', true)
end
if doplot.VoltageMesh_OHC_rel
    plot_voltage_mesh('OHC', true)
    plot_voltage_mesh('OHC_MET', true)
end
if doplot.VoltageMesh_IHC
    plot_voltage_mesh('IHC', false)
    plot_voltage_mesh('IHC_MET', false)
end
if doplot.VoltageMesh_OHC
    plot_voltage_mesh('OHC', false)
    plot_voltage_mesh('OHC_MET', false)
end


%% Channels HC

    function plot_channels_mesh(channel, ss)
        
        V_ss = subsampleVoltage(fChannels, channel.name, 1, mnaopt, plotopt, 'channels');
        V_t = subsampleVoltage(fChannels, channel.name, ts, mnaopt, plotopt, 'channels');

        if ss == true
            V_t = V_t - V_ss;
            V0 = 0;
        else
            V0 = [];
        end

        hfig = mna.plot.voltage_mesh_one(Xgrid, Tgrid, V_t, t0, tf, plotopt, V0);
        
        zlabel('p Open');
        
        title(channel.name)
        
        % save figure        
        if runopt.save_figures
            fname = sprintf('Channels_%s', channel.name);
            plotopt.save_figure(hfig, runopt.figureSaveDir, fname);
        end
    end


if doplot.channels_mesh
    channels = mnaopt.channels;
    for i = 1:numel(channels)
        plot_channels_mesh(channels(i), false)
    end
end
if doplot.channels_mesh_rel
    channels = mnaopt.channels;
    for i = 1:numel(channels)
        plot_channels_mesh(channels(i), true)
    end
end

%% Mech

    function plot_mech_mesh(v, ss)
        
        V_ss = subsampleVoltage(fMech, v, 1, mechopt, plotopt, 'mech');
        V_t = subsampleVoltage(fMech, v, ts, mechopt, plotopt, 'mech');

        if ss == true
            V_t = V_t - V_ss;
            V0 = 0;
        else
            V0 = [];
        end

        hfig = mna.plot.voltage_mesh_one(Xgrid_mech, Tgrid_mech, V_t, t0, tf, plotopt, V0);
        
        zlabel('d [nm]');
        
        title(v)
        
        % save figure        
        if runopt.save_figures
            fname = sprintf('%s', v);
            plotopt.save_figure(hfig, runopt.figureSaveDir, fname);
        end
    end

vars = {'BMx', 'TMx', 'BMv', 'TMv'};
for i = 1:numel(vars)
    vv = vars{i};
    if doplot.(['mech_', vv])
        plot_mech_mesh(vv, false)
    end
    
    if doplot.(['mech_profile_', vv])
        mna.plot.profile_one_graph(xx, ...
            @(i) fMech.(vv), ...
            {vv}, mnaopt, plotopt, runopt, t1);
        title(vv);
        ylabel('displ');
    end
    
end

%% AMPL

function plot_ampl_mesh(v)

    V_ss = subsampleVoltage(fChannels, 'OHC_fast', 1, mnaopt, plotopt, 'channels');
    V_t = subsampleVoltage(fChannels, 'OHC_fast', ts, mnaopt, plotopt, 'channels');

    V_t_channel = V_t - V_ss;
    
    V_ss = subsampleVoltage(fVolt, 'OHC', 1, mnaopt, plotopt, 'voltage');
    V_t = subsampleVoltage(fVolt, 'OHC', ts, mnaopt, plotopt, 'voltage');
    
    V_t_ohc = V_t - V_ss;
    
    V_ss = subsampleVoltage(fMech, v, 1, mechopt, plotopt, 'mech');
    V_t = subsampleVoltage(fMech, v, ts, mechopt, plotopt, 'mech');

    V_t_mech = V_t - V_ss;
    
    V0 = 0;
    
    % VAR = V_t_mech .* V_t_channel;
    VAR = V_t_mech .* V_t_ohc;
    
    hfig = mna.plot.voltage_mesh_one(Xgrid_mech, Tgrid_mech, VAR, t0, tf, plotopt, V0);

    zlabel('ampl');

    title(v)

    % save figure        
    if runopt.save_figures
        fname = sprintf('%s', v);
        plotopt.save_figure(hfig, runopt.figureSaveDir, fname);
    end
end

if false
    plot_ampl_mesh('BMv')
end
    


%% Voltage

if doplot.VoltageMesh
    
    names = {'IHC_ic', 'OHC_ic', 'ScalaMedia', 'IHC_ec', 'OHC_ec', 'SpiralLigament'};
    fullnames = {
        'IHC intracellular', ...
        'OHC intracellular', ...
        'Scala Media', ...
        'IHC extracellular', ...
        'OHC extracellular', ...
        'Spiral Ligament'};
    vars = struct();
    
    for i = 1:numel(names)
        vars.(names{i}) = struct( ...
            'name', fullnames{i}, ...
            'data', subsampleVoltage(fVolt, names{i}, ts, mnaopt, plotopt, 'voltage'), ...
            'steadystate', subsampleVoltage(fVolt, names{i}, 1, mnaopt, plotopt, 'voltage'));
    end    
    
    mna.plot.voltage_mesh( Xgrid, Tgrid, vars, t0, tf, plotopt, runopt )
end

if doplot.CochlearMicrophonic
    
    V_ss = subsampleVoltage(fVolt, 'N1', 1, mnaopt, plotopt, 'voltage');
    V = subsampleVoltage(fVolt, 'N1', ts, mnaopt, plotopt, 'voltage');
    V = V - V_ss;
    
    BMlength = 28.8;
    BMatt = 13.38;
    BM1el = 2.68;
    
    el_width = 0.4; % mm
    
    BM18el = BM1el + BMatt;    
    BM22el = BM18el + 2.5;
    
    el118 = linspace(BM1el, BM18el, 18-1+1);
    el1822 = linspace(BM18el, BM22el, 22-18+1);
    
    electrodes = [el118(1:end-1), el1822] / BMlength;
    
    %offset = 2.0 /BMlength;
    offset = 0.01 /BMlength;
    
    electrodes = electrodes + offset;
       
    zz = flip(electrodes(2:2:end));
    
    for kkk = 1:2
        figure
        hold on

        surf(Xgrid, Tgrid, V, ...
        plotopt.surf_options{:});

        xlabel('BM position');
        ylabel('Time [ms]');
        ylim([t0.ms, tf.ms]);
        zlabel('Voltage [mV]');
        title('node 1 - CM');

        M = max(abs(V(:)));
        caxis([-M, M])

        colormap(plotopt.centered_colormap)
        view(plotopt.view);

        set(gca, 'color', 0.95*ones(1,3));

        for z = zz
            plot3([z,z],[0,tf.ms],[0,0],'k-')
        end
        
        if kkk == 1
            set(gca, 'view', [0,90])
        end
        if kkk == 2
            set(gca, 'view', [-8,16])
        end
    end
    
    VVV = max(V(:)) - min(V(:));
    
    figure
    hold on
    
    ii = zeros(size(zz));
    kk = 22;
    for j = 1:numel(zz)
        z = zz(j);

        % only for uniform sampling    
        % i = round(z*mnaopt.Numstacks);    
        % di = el_width/2/BMlength * mnaopt.Numstacks;    
        % di = round(di);
        % i0 = i - di;
        % i1 = i + di;

        i = find(mnaopt.xgrid >= z, 1);

        i0 = find(mnaopt.xgrid >= z - el_width/2/BMlength, 1);
        i1 = find(mnaopt.xgrid >= z + el_width/2/BMlength, 1);

        ii(j) = i;

        VV = mean(V(:,i0:i1), 2) / VVV;
        
        plot(T, kk + 6*VV);        
        kk = kk - 2;
    end
    
    yticks([2:2:22])
    
    ylabel('node 1 - CM amplitude by electrode')
    xlabel('time (ms)')
    
    figure
    hold on
    
    VV = (max(V, [], 1) - min(V, [], 1)) / VVV;
    
    plot(xx, VV);
    plot(xx(ii), VV(:,ii), 'ko');
    
    xlabel('rel. distance from base')
    ylabel('rel. magnitude')
    title('CM peak-to-peak')
    
    legend({'node 1 voltage','CI electrodes'})
    
    figure
    hold on
    
    f = stimulus.eval_fcn('ms');
    S = f(ts);
    S = S(:);
        
    phase_lag_hilbert = zeros(size(xx));
    for i = 1:numel(xx)
        VAR = V(:,i);
        VAR = VAR(:);
        
        phase_lag_hilbert(i) = mean(unwrap(angle(hilbert(VAR))) - unwrap(angle(hilbert(S))));
    end
    
    ind = xx < 0.75;
    
    ph = unwrap(phase_lag_hilbert(ind));
    ph = 180/pi*ph;
    ph = ph - ph(1);
    
    plot(xx(ind), ph);
    plot(xx(ii), ph(ii), 'ko');
    
    xlim([0,1])
    
    % xBM = resFiles.fBM.x;
    xBM = mnaopt.xgrid;
    
    BMx = resFiles.fMech.BMx(ts,:);
    
    phase_lag_hilbert = zeros(size(xx));
    for i = 1:numel(xx)
        
        [~, xind] = min(abs(xx(i) - xBM));
        S = BMx(:,xind);
        
        S = S(:);
        
        VAR = V(:,i);
        VAR = VAR(:);
        
        phase_lag_hilbert(i) = mean(unwrap(angle(hilbert(VAR))) - unwrap(angle(hilbert(S))));
    end
    
    ind = xx < 0.75;
    
    ph = unwrap(phase_lag_hilbert(ind));
    ph = 180/pi*ph;
    ph = ph - ph(1);
    
    plot(xx(ind), ph);
    plot(xx(ii), ph(ii), 'ko');
    
    xlabel('rel. distance from base')
    ylabel('phase')
    title('CM phase')
    
    xlim([0,1])
    
    legend({'acoustic phase diff.','CI electrodes','BM displ. phase diff.','CI electrodes'})
    
    ax = gca;
    
    x0 = BM18el / BMlength;
    x1 = BM22el / BMlength;
    
    a = fill([x0 x0 x1 x1], [ax.YLim(2) ax.YLim(1) ax.YLim(1) ax.YLim(2)], [0,0,0]);
    a.FaceAlpha = 0.1;
    a.EdgeAlpha = 0.0;
    
    
end

%% Voltage Steady State Hair Cell

if doplot.VoltageSteadyState_HC
    mna.plot.voltage_profile_hc(xx, fVolt, mnaopt, plotopt, runopt,t1);
end

%% Voltage Steady State

if doplot.VoltageSteadyState
    mna.plot.voltage_profile(xx, fVolt, mnaopt, plotopt, runopt,t1);
end

%% Voltage Profile

% =========================================================================
% RESISTORS
ele = mnaopt.circuit.select('type resistor and special radial');

% -------------------------------------------------------------------------
if doplot.voltage_profile_one_graph_resistors
    
    resistors = {'IHC', 'OHC'};
    
    hfig = mna.plot.profile_one_graph(xx, ...
        @(i) sim_to_mV * getVoltage(fVolt, resistors{i}, mnaopt, 'transformed'), ...
        resistors, mnaopt, plotopt, runopt, t1);
    title('hair cell potentials');
    ylabel('Voltage [mV]');
    
    if runopt.save_figures
        mySaveFig(hfig, "voltage_profile_one_graph_HC", [], runopt, ...
                    'relativeDataPath', 'img/MNA/');
    end
    
    resistors = {ele.name};
    % resistors = {'RAI', 'RBI', 'RSL', 'RStV', 'RAO', 'RBO', 'ROC'};
    
    resistors = [resistors, setdiff(resistors, {ele.name})];
    
    hfig = mna.plot.profile_one_graph(xx, ...
        @(i) sim_to_mV * getVoltage(fVolt, resistors{i}, mnaopt, 'transformed'), ...
        resistors, mnaopt, plotopt, runopt, t1);
    title('resistors');
    ylabel('Voltage [mV]');
    
    if runopt.save_figures
        mySaveFig(hfig, "voltage_profile_one_graph_resistors", [], runopt, ...
                    'relativeDataPath', 'img/MNA/');
    end
    
%     resistors = {'RRM', 'RSVSL', 'RBM'};
%     
%     hfig = mna.plot.profile_one_graph(xx, ...
%         @(i) sim_to_mV * getVoltage(fVolt, resistors{i}, mnaopt, 'transformed'), ...
%         resistors, mnaopt, plotopt, runopt, t1);
%     title('resistors');
%     ylabel('Voltage [mV]');
    
    

end
% -------------------------------------------------------------------------
if doplot.voltage_profile_subplot_resistors
    
    resistors = {'IHC', 'OHC', ele.name};
    
    mna.plot.voltage_profile_subplot(xx, ...
        @(i) sim_to_mV * getVoltage(fVolt, resistors{i}, mnaopt, 'transformed'), ...
        resistors, mnaopt, plotopt, runopt, t1);
end


% =========================================================================
% LONGITUDINAL RESISTORS

% steady state value of BM and OHC stereocilia displacement
BM_displ = zeros(mnaopt.Numstacks,1);
OHC_cilia = zeros(mnaopt.Numstacks,1);

% resting potentials of the inner and outer hair cells before iterative setting
vihc = mnaopt.IHC_V_rest(mnaopt.xgrid(:));
vohc = mnaopt.OHC_V_rest(mnaopt.xgrid(:));

dependencies = { ...
    'xpos', mnaopt.xgrid(1:end-1), ...
    'BM_displ', BM_displ, ...
    'OHC_cilia', OHC_cilia, ... 
    'vihc', vihc, ...
    'vohc', vohc, ...
    'vihc_ss', vihc, ...
    'vohc_ss', vohc};

ele = mnaopt.circuit.select('type resistor and special longitudinal');
long_resistors = {ele.name};

for i = 1:numel(long_resistors)
    G{i} = mnaopt.circuit.get_element_by_name(long_resistors{i}).value(dependencies{:});
end

% -------------------------------------------------------------------------
if doplot.voltage_profile_one_graph_long_resistors
    mna.plot.profile_one_graph(xx(1:end-1), ...
        @(i) sim_to_mV * getVoltage(fVolt, long_resistors{i}, mnaopt, 'transformed'), ...
        long_resistors, mnaopt, plotopt, runopt, t1)
    title('longitudinal resistors');
    ylabel('Voltage [mV]');
end
% -------------------------------------------------------------------------
if doplot.voltage_profile_subplot_long_resistors
    mna.plot.voltage_profile_subplot(xx(1:end-1), ...
        @(i) sim_to_mV * getVoltage(fVolt, long_resistors{i}, mnaopt, 'transformed'), ...
        long_resistors, mnaopt, plotopt, runopt,t1);
end
% -------------------------------------------------------------------------
if doplot.current_profile_one_graph_long_resistors
    mna.plot.profile_one_graph(xx(1:end-1), ...
        @(i) sim_to_mV * G{i} .* getVoltage(fVolt, long_resistors{i}, mnaopt, 'transformed'), ...
        long_resistors, mnaopt, plotopt, runopt,t1)
    title('longitudinal resistors');
    ylabel('Current [mA]');
end


% =========================================================================
% NODES
nodes = cellfun(@(x) sprintf('N%d', x), num2cell(1:mnaopt.circuit.num_node), 'UniformOutput', false);

% -------------------------------------------------------------------------
if doplot.voltage_profile_one_graph_nodes
    mna.plot.profile_one_graph(xx, ...
        @(i) sim_to_mV * getVoltage(fVolt, nodes{i}, mnaopt, 'transformed'), ...
        nodes, mnaopt, plotopt, runopt, t1)
    title('nodes');
    ylabel('Voltage [mV]');
end
% -------------------------------------------------------------------------
if doplot.voltage_profile_subplot_nodes
    mna.plot.voltage_profile_subplot(xx, ...
        @(i) sim_to_mV * getVoltage(fVolt, nodes{i}, mnaopt, 'transformed'), ...
        nodes, mnaopt, plotopt, runopt, t1);
end

%% Show figures

if ~isempty(plotopt.figure_handles)
    arrayfun(@(h) set(h, 'Visible', 'On'), plotopt.figure_handles);
end

%% Subsample voltage


    function V = subsampleVoltage(fVolt, var, ts, mnaopt, plotopt, sw)
        switch sw
            case 'voltage'
                [ ~, sel, name, operand, fun ] = getVoltage( [], var, mnaopt, 'transformed');
            case 'channels'
                name = var;
                operand = [];
                sel = 1:mnaopt.Numstacks;
            case 'mech'
                name = var;
                operand = [];
                sel = 1:mnaopt.Numstacks;
        end
        
        if isempty(operand)
            V = plotopt.subsample(fVolt, name, ts, sel);
        else
            assert(numel(name) == numel(sel));
            for k = numel(name) : -1 : 1
                v{k} = plotopt.subsample(fVolt, name{k}, ts, sel{k});
            end
            V = fun(v{:});
        end
        
        switch sw
            case 'voltage'
                V = V * sim_to_mV;
            case 'channels'
                % V = V;
        end
                    
    end

end