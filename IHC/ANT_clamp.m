function [ ANTStatistics ] = ANT_clamp( stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt, plotopt, args )
arguments
    stimulus
    topt
    midopt
    mechopt
    mnaopt
    antopt
    hhopt
    runopt
    opt
    memopt
    paropt
    plotopt

    args.init_flag = false
    args.xpos = []
end
%ANT
disp_title( 'Auditory Nerve Transduction Simulation', '*' );

%% Init
disp_title( 'Synapse INIT' );

% [potential, displacement, TimeVoltage, Signal, TimeSignal, antopt, runopt] = ...
%     iniANT( stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt );

% xpos = antopt.positionsToExcite;

%%

topt.compute = copy(topt.total);

[cx, cf] = characteristic_frequency_model('Li');

if isempty(args.xpos)
    xpos = cx(stimulus.frequency.Hz);
else
    xpos = args.xpos;
end

antopt.slicesToExcite = 1;
antopt.positionsToExcite = xpos;

fs = stimulus.samplingFrequency;

T0 = stimulus.zeroDuration.s;
T2 = stimulus.fadeDuration.s;

T1 = stimulus.Duration.s - T0 - T2;

[t, ~] = single_IHC.step_signal([T0, T1, T2], [0, 0, 0], "fs", fs.Hz);

displacement = zeros(size(t));

%%

channels = mnaopt.IHC_basolateral_channels;

for i = 1:numel(channels)
    G_grad = channels(i).parameters.x_grad.G;
    EK_grad = channels(i).parameters.x_grad.EK;
    channels(i).parameters.G = linf(G_grad(1), G_grad(2), xpos);
end

out = single_IHC.main(xpos, displacement * 1e-9, fs.Hz, ...
    "EP", 90e-3, ...
    "EK", linf(EK_grad(1), EK_grad(2), xpos), ...
    "channels", channels);

Vm = out.Vm;


%%
V0 = mean(Vm);

Vs = stimulus.amplitude*1e-3;

V1 = V0 + Vs;

[t, step] = single_IHC.step_signal([T0, T1, T2], [V0, V1, V0], "fs", fs.Hz);
step = movingAverageFilter(0.1e-3*fs, step);

t = t(:);
step = step(:);

Signal = zeros(size(t));
TimeSignal = Time(t, 's');
TimeVoltage = Time(t, 's');

%% Preprocessing
% * interpolate the potential (Voltage) Time and Signal
%   to respect the required sampling
% * scale voltage by provided factor

disp('Using rand modification of IHC potential')
potential = step *1e3;

disp_title( 'Preprocessing' );

[potential, time, Signal] = ANTPreprocess( potential, TimeVoltage, Signal, TimeSignal, antopt, topt );

%% Synapse IC

use_IC = false;

if use_IC && ~args.init_flag
    extra_dirs = {runopt.read_only_cache};
            
    [runopt, ~] = findResultExtended(extra_dirs, 'synapse_ic', ...
            [], [], [], mechopt, mnaopt, antopt, [], runopt, opt);
    
    cache_ic_file = fullfile(runopt.path.synapse_ic, 'synapse_ic.mat');
    if ~exist(runopt.path.synapse_ic, 'dir')
        mkdir(runopt.path.synapse_ic)
    end
    
    if ~runopt.recalculate.synapse_ic ...
            && runopt.found.synapse_ic ...
            && exist(cache_ic_file, 'file')
        dbprintf('Loading initial conditions from a file:\n -- %s.\n', cache_ic_file);
        IC = load(cache_ic_file);
    elseif runopt.no_create_only_load_synapse_ic
        
        error('We require initial conditions to exist.')
        
    else
    
        [topt_ic, stimulus_ic] = devopts.stimulus( ...
                antopt.samplingFrequency, Time(5000, 'ms'), ...
                'special', 'zero');

        antopt_ic = copy(antopt);
        antopt_ic.numberOfRepetitions = 1;

        runopt_ic = copy(runopt);

        runopt_ic.do.nerve = false;
        runopt_ic.plot.synapse = false;
        runopt_ic.plot.nerve = false;

        runopt_ic.analysisStart = Time(100, 'ms');
    
        ANTStatistics_IC = ANT_single(stimulus_ic, topt_ic, ...
            copy(midopt), copy(mechopt), copy(mnaopt), antopt_ic, copy(hhopt), ...
            runopt_ic, copy(opt), copy(memopt), copy(paropt), copy(plotopt), ...
            'xpos', xpos, ...
            'init_flag', true);
        f = load(ANTStatistics_IC.main.Properties.Source);
        
        tmp = [f.FreePool.TimeAverage];
        q0 = round([tmp.MEAN]);
        
        tmp = [f.CaConc.TimeAverage];
        C0 = [tmp.MEAN];

        tmp = [f.ReStore.TimeAverage];
        w0 = round([tmp.MEAN]);

        m0 = f.deterministic_stats.m.statistics.tspan.mean;
        
        IC = struct( ...
            'q0', q0, ...
            'C0', C0, ...
            'w0', w0, ...
            'm0', m0);

        save(cache_ic_file, '-struct', 'IC', '-v7.3');

        % Create SUCCESS indicator
        fid = fopen(fullfile(runopt.path.synapse_ic, 'SUCCESS'),'w');
        fclose(fid);
    end
else
    % IC for calculating IC
    IC = [];
end

antopt.initialConditions = IC;

%% Synapse IC
disp_title( 'Synapse Simulation' );

% get synapse result
[SynResFile, runopt, antopt] = handleSynapseResult( potential, time, stimulus, topt, midopt, mechopt, mnaopt, antopt, runopt, paropt, opt, memopt );

% save settings
% runopt = saveANTSettings( hhopt, runopt );


%% Nerve
if runopt.do.nerve
    disp_title( 'Nerve Simulation' );

    % get nerve result
    [NerveResFile, runopt, hhopt] = handleNerveResult( SynResFile, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, paropt, opt, memopt );
else
    NerveResFile = [];
end

%% postprocess result

postprocess_windows = {'full', 'main'};

if runopt.do.ant_postprocess
    for i = 1:numel(postprocess_windows)
        w = postprocess_windows{i};
        StatFileName = fullfile(runopt.path.synapse, ...
            sprintf('ANTStatistics_%s.mat', w));
        
        doANTPostprocess = false;
        
        if runopt.recalculate.ant_postprocess || runopt.recalculate.synapse
            doANTPostprocess = true;
        end
    
        if ~exist(StatFileName,'file')
            doANTPostprocess = true;
        end
    
        if exist(StatFileName,'file')
            load(StatFileName, 'num_replicas_synapse', 'num_replicas_nerve');
            if ~exist('num_replicas_synapse', 'var') || ~exist('num_replicas_nerve', 'var')
                warning('Variable num_replicas_synapse or num_replicas_nerve not found in ANTStatistics.mat, recalculating');
                doANTPostprocess = true;
            else
                if num_replicas_synapse < antopt.numberOfRepetitions
                    doANTPostprocess = true;
                end
                if num_replicas_nerve < antopt.numberOfRepetitions
                    doANTPostprocess = true;
                end
            end
        end
        
        if doANTPostprocess
            ANTStatistics.(w) = ANTPostprocess( ...
                SynResFile, NerveResFile, potential, time, Signal, ...
                antopt, hhopt, runopt, stimulus, memopt, ...
                'time_span', w);
        else
            ANTStatistics.(w) = matfile(StatFileName,'Writable', true);
        end
    end
else
    ANTStatistics = [];
end

%% Plotting

% postprocess_windows_to_plot = {'full', 'main'};
postprocess_windows_to_plot = {'full'};

if runopt.plot.ant_postprocess
    if isempty(ANTStatistics)
        warning('Can''t plot ANT Postprocess -- variable empty. Make sure runopt.do.ant_postprocess is set to true.')
    else
        for i = 1:numel(postprocess_windows_to_plot)
            w = postprocess_windows_to_plot{i};
            ANTPostprocessPlotting( ANTStatistics.(w), SynResFile, NerveResFile, potential, Time, Signal, antopt, hhopt, runopt, plotopt, stimulus, 1 )
        end
    end
end

if runopt.plot.ant || runopt.plot.synapse || runopt.plot.nerve
    
    if ischar(antopt.plotSynapseSimulation)
        if strcmp(antopt.plotSynapseSimulation, 'all')
            L = 1 : antopt.numberOfRepetitions;
        end
    elseif isnumeric(antopt.plotSynapseSimulation)
        L = antopt.plotSynapseSimulation;
    else
        L = [];
    end
       
    SEL = [];
    for i = L
        if any(strcmp(SynResFile{i}.SR, runopt.SRsToPlot)) || any(strcmp('all', runopt.SRsToPlot))
            SEL(end+1) = i;
        end
    end
end

% plot ANT
if runopt.plot.ant 
    if numel(antopt.slicesToExcite) > 1
        ANTplotting( SynResFile(SEL), NerveResFile(SEL) )
    end
end

% plot synapse results
if runopt.plot.synapse
    maxslice_relative = slicesToExcite({'max'}, mechopt.Numstacks, SynResFile{1}.V, [], mnaopt.xgrid);

    fprintf('Plotting slice %d ...\n', antopt.slicesToExcite(maxslice_relative))
    for i = 1:numel(postprocess_windows_to_plot)
        w = postprocess_windows_to_plot{i};
        SynapsePlotting_avg( potential, time, ANTStatistics.(w), SynResFile(SEL), NerveResFile(SEL), antopt, hhopt, opt, runopt, plotopt, maxslice_relative );
    end

%     disp_title( 'Synapse Plotting' );
%     if numel(SEL) > 5
%         SEL = SEL(1:5);
%         warning('Attempting to plot too many simulations, the number reduced to 5')
%     end
%     
%     fprintf('Plotting slice %d ...\n', antopt.slicesToExcite(maxslice_relative))
%     for i = 1:numel(postprocess_windows)
%         w = postprocess_windows{i};
%         SynapsePlotting( potential, time, ANTStatistics.(w), SynResFile(SEL), NerveResFile(SEL), antopt, hhopt, opt, runopt, plotopt, maxslice_relative );
%     end
end

% plot nerve results
if runopt.plot.nerve
    disp_title( 'Nerve Plotting' );
    if numel(SEL) > 5
        SEL = SEL(1:5);
        warning('Attempting to plot too many simulations, the number reduced to 5')
    end
    for i = 1:numel(SEL)
        nerveplotting(NerveResFile{SEL(i)}, stimulus, topt, mechopt, mnaopt, antopt, hhopt, runopt, opt, plotopt, memopt);
    end
end

end
