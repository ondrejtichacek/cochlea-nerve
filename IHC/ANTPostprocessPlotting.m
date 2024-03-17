function ANTPostprocessPlotting( src, SynResFile, NerveResFile, Voltage, time, Signal, antopt, hhopt, runopt, plotopt, stimulus, i )
%ANTPOSTPROCESSPLOTTING

%% INIT

% Currenlty we don't calculate statistics on N_Na and Istim - needs to be
% reimplemented.

% [Cleft, FreePool, ReStore, NTRate, Volt, N_Na, Istim, CaConc] = loadv(src, ...
    % 'Cleft', 'FreePool', 'ReStore', 'NTRate', 'Volt', 'N_Na', 'Istim', 'CaConc');

[Cleft, FreePool, ReStore, NTRate, Volt, CaConc] = loadv(src, ...
    'Cleft', 'FreePool', 'ReStore', 'NTRate', 'Volt', 'CaConc');

[Signal_Toolbox_licenceAvailable, ~] = license('checkout','Signal_Toolbox');

if ~Signal_Toolbox_licenceAvailable
    warning('Signal toolbox not available, skipping some parts of analysis!');
end

SR = getActiveSR(antopt);

hfig = [];

%% Restrict to i-th cross-section

warning('Plotting %d-th relative crossection');

Cleft = Cleft(i);
FreePool = FreePool(i);
ReStore = ReStore(i);
NTRate = NTRate(i);
Volt = Volt(i);
% N_Na = N_Na(i);
% Istim = Istim(i);
CaConc = CaConc(i);

%% Restrict to analysis time frame

analysisStartTime = runopt.analysisStart;
analysisStartTimeIndex = min([length(time), find(analysisStartTime <= time,1)]);
warning('Analysis domain striped to < %s , %s >.', time(analysisStartTimeIndex), time(end));

% tmp fix
analysisStartTimeIndex = analysisStartTimeIndex + 1;

time = time(analysisStartTimeIndex:end);
Voltage = Voltage(analysisStartTimeIndex:end,i);


%% Hilbert Transform plots

hfig(end+1) = plotopt.figure;

FACTOR = 100;
warning('HT plots subsampled by factor %d', FACTOR);

% -------------------------------------------------------------------------
subplot(2,2,1)
var = src.SignalAmplitude(analysisStartTimeIndex:end, :);
[S, T] = intelligentSubsampling( ...
    var, time.ms, FACTOR, ...
    'Npeaks', round(length(time)/FACTOR), ...
    'Threshold', 10);

plot(T, S);

xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')

set(gca, ...
    'XLim', [0,time(end).ms]);

% -------------------------------------------------------------------------
subplot(2,2,2)
var = src.InstantaneousFrequency(analysisStartTimeIndex:end, :);
[S, T] = intelligentSubsampling( ...
    var, time(1:end-1).ms, FACTOR, ...
    'Npeaks', round(length(time)/FACTOR), ...
    'Threshold', 5);

plot(T, S);

xlabel('Time (ms)')
ylabel('Frequency (Hz)')

set(gca, ...
    'XLim', [0,time(end).ms], ...
    'YLim', [0,2*max([stimulus.frequency.Hz])]);

% -------------------------------------------------------------------------
subplot(2,2,3:4)
var = src.SignalPhase(analysisStartTimeIndex:end, :);
[S, T] = intelligentSubsampling( ...
    var, time.ms, FACTOR);

plot(T*1e3, S);

xlabel('Time (ms)')
ylabel('Phase (rad)')

set(gca, ...
    'XLim', [0,100], ...
    'YLim', [-pi,pi], ...
    'YTick', [-pi: pi/2 : pi], ...
    'YTickLabel', {'$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'});

% -------------------------------------------------------------------------
if runopt.save_figures
    mySaveFig( hfig(end), ...
        'SignalHilberTransform', ...
        runopt.path.synapse, ...
        runopt.path.synapse, ...
        'width', '0.9\textwidth', ...
        'height', '0.5\textwidth' );
end
    
%% Fourier transform plots

% Voltage
[fV, VV, VVn] = plotSignalFT( ...
    Voltage - mean(Voltage), ...
    hhopt.samplingFrequency.Hz );

% Signal
[fS, SS, SSn] = plotSignalFT( ...
    Signal - mean(Signal), ...
    hhopt.samplingFrequency.Hz );

[~,ind] = max(SS);

MainSignalFrequency = fS(ind)

FACTOR = 10;
warning('FT plots subsampled by factor %d', FACTOR);
subsamplingFT = unique([1:FACTOR:length(fV),length(fV)]);

for i = 1:length(SR)
    if any(strcmp(SR{i}, runopt.SRsToPlot)) || any(strcmp('all', runopt.SRsToPlot))
        
        % =================================================================
        % Cleft
        for k = find(strcmp({Cleft.sr}, SR{i}))
            [f,YC,YCn] = plotSignalFT( ...
                Cleft(k).MEAN - mean(Cleft(k).MEAN), ...
                hhopt.samplingFrequency.Hz);

            LegendText = {sprintf('FT of Mean Cleft %s',SR{i}), 'FT of IHC Voltage'};

            hfig(end+1) = myFTplot(f, YC, fV, VV, ...
                LegendText, sprintf('== %s ==',SR{i}), ...
                subsamplingFT, plotopt);

            if runopt.save_figures
                mySaveFig( hfig(end), ...
                    sprintf('SynMeanCleftFT_%s.fig',SR{i}), ...
                    runopt.path.synapse, ...
                    runopt.path.synapse, 'width', '0.7\textwidth' );
            end
        end

        % =================================================================
        % Volt
        for k = find(strcmp({Volt.sr}, SR{i}))
            if runopt.do.nerve
                [f,YY,YYn] = plotSignalFT( ...
                    Volt(k).MEAN - mean(Volt(k).MEAN), ...
                    hhopt.samplingFrequency.Hz );

                LegendText = {sprintf('FT of Mean Nerve Volt %s',SR{i}), 'FT of IHC Voltage'};

                hfig(end+1) = myFTplot(f, YY, fV, VV, ...
                    LegendText, sprintf('== %s ==',SR{i}), ...
                    subsamplingFT, plotopt);

                if runopt.save_figures
                    mySaveFig( hfig(end), ...
                        sprintf('NerveMeanVoltFT_%s',SR{i}), ...
                        runopt.path.nerve, ...
                        runopt.path.nerve, 'width', '0.7\textwidth' );
                end
            end
        end
        % =================================================================
        % Joined    
        if runopt.do.nerve
            LegendText = {  sprintf('Mean Synapse Cleft %s',SR{i}), ...
                            sprintf('Mean Nerve Volt %s',SR{i}), ...
                            'IHC Voltage', ...
                            'Accoustic Input Signal'};

            hfig(end+1) = myJoinedFTplot( {f(subsamplingFT),f(subsamplingFT),f(subsamplingFT),f(subsamplingFT)}, ...
                             {YCn(subsamplingFT), YYn(subsamplingFT)}, ...
                             {VVn(subsamplingFT), SSn(subsamplingFT)}, ...
                             LegendText, plotopt);

             if runopt.save_figures
                 mySaveFig( hfig(end), ...
                    sprintf('joinedFT_%s',SR{i}), ...
                    runopt.path.nerve, ...
                    runopt.path.nerve, 'width', '0.9\textwidth', 'height', '0.3\textwidth' );
             end
        end
    end
end
    
%% Centered yy figures

FACTOR = 50;
warning('Other plots subsampled by factor %d', FACTOR);
subsampling = unique([1:FACTOR:length(time),length(time)]);

for i = 1:length(SR)
    if any(strcmp(SR{i}, runopt.SRsToPlot)) || any(strcmp('all', runopt.SRsToPlot))
        
        % =================================================================
        % Cleft
        for k = find(strcmp({Cleft.sr}, SR{i}))
            hfig(end+1) = plotopt.figure;
            myCenteredPlot( ...
                time(subsampling).ms, Cleft(k).MEAN(subsampling), ...
                time(subsampling).ms, Voltage(subsampling), ...
                0, Voltage(1), ...
                'Time (ms)', 'Mean NT Vesicles (1)', 'IHC Voltage (mV)');       
            
            if runopt.save_figures
                mySaveFig( hfig(end), ...
                    sprintf('SynMeanCleft_%s',SR{i}), ...
                    runopt.path.synapse, ...
                    runopt.path.synapse, 'width', '0.9\textwidth', 'height', '0.3\textwidth' );
            end
        end
        
        % Istim
        % for k = find(strcmp({Istim.sr}, SR{i}))
        %     hfig(end+1) = plotopt.figure;
        %     myCenteredPlot( ...
        %         time(subsampling).ms, Istim(k).MEAN(subsampling), ...
        %         time(subsampling).ms, Voltage(subsampling), ...
        %         0, Voltage(1), ...
        %         'Time (ms)', 'Mean Nerve Istim (mA)', 'IHC Voltage (mV)');       
        % 
        %     if runopt.save_figures
        %         mySaveFig( hfig(end), ...
        %             sprintf('NerveMeanIstim_%s',SR{i}), ...
        %             runopt.path.nerve, ...
        %             runopt.path.nerve, 'width', '0.9\textwidth', 'height', '0.25\textwidth' );
        %     end
        % end
        
        % Volt
        for k = find(strcmp({Volt.sr}, SR{i}))
            hfig(end+1) = plotopt.figure;
            myCenteredPlot( ...
                time(subsampling).ms, Volt(k).MEAN(subsampling), ...
                time(subsampling).ms, Voltage(subsampling), ...
                0, Voltage(1), ...
                'Time (ms)', 'Mean Nerve Volt (mV)', 'IHC Voltage (mV)');       

            if runopt.save_figures
                mySaveFig( hfig(end), ...
                    sprintf('NerveMeanVolt_%s',SR{i}), ...
                    runopt.path.nerve, ...
                    runopt.path.nerve, 'width', '0.9\textwidth', 'height', '0.25\textwidth' );
            end
        end
    
        % =================================================================
        % NT Rate figure
        for k = find(strcmp({NTRate.sr}, SR{i}))
            hfig(end+1) = plotopt.figure;
            plot(time(subsampling).ms, NTRate(k).MEAN(subsampling));
            % title('NT Rate')
            xlabel('Time (ms)')
            ylabel('NT Rate (Hz)')

            if runopt.save_figures
                mySaveFig( hfig(end), ...
                    sprintf('SynMeanNTRate_%s',SR{i}), ...
                    runopt.path.synapse, ...
                    runopt.path.synapse, 'width', '0.9\textwidth', 'height', '0.2\textwidth' );
            end
        end

        % =================================================================
        % Stores figure
        k1 = find(strcmp({FreePool.sr}, SR{i}));
        k2 = find(strcmp({ReStore.sr}, SR{i}));        
        assert(numel(k1) == numel(k2));
        
        for k = 1:numel(k1)
            hfig(end+1) = plotopt.figure;
            hold on
        
%         for j = 1:5
%             f = SynResFile{j};
%             if strcmp(SR{i}, f.SR)
%                 
%                 q = f.q;
%                 
%                 FACTOR = 1000;
%                 
%                 [S, T] = intelligentSubsampling( ...
%                     q, Time*1e3, FACTOR, ...
%                     'Npeaks', round(length(Time)/FACTOR), ...
%                     'Threshold', 10);
%                 
%                 plot(T, S, 'Color', 0.9*[1,1,1]);
%             end
%         end
        
            plot(time(subsampling)*1e3, FreePool(k1(k)).MEAN(subsampling));
            plot(time(subsampling)*1e3, ReStore(k2(k)).MEAN(subsampling));

            % title('NT dynamics')
            xlabel('Time (ms)')
            ylabel('NT Vesicles (1)')

            legend({ sprintf('Mean q (%s)', SR{i}), sprintf(' Mean w (%s)', SR{i}) });
            
            if runopt.save_figures
                mySaveFig( hfig(end), ...
                    sprintf('SynMeanStores_%s',SR{i}), ...
                    runopt.path.synapse, ...
                    runopt.path.synapse, 'width', '0.9\textwidth', 'height', '0.35\textwidth' );
            end
        end
    
        % =================================================================
        % InterTimes Histogram
        k1 = find(strcmp({Cleft.sr}, SR{i}));
        k2 = find(strcmp({Volt.sr}, SR{i}));        
        assert(numel(k1) == numel(k2));
        
        for k = 1:numel(k1)
    
            hfig(end+1) = plotopt.figure;
            hold on

            histogram( cat(1, Cleft(k1(k)).InterTimes{:}), 'Normalization', 'pdf' );
            histogram( cat(1, Volt(k2(k)).InterTimes{:}), 'Normalization', 'pdf' );

            title('Inter-peak-time Histogram')
            xlabel('Time (ms)')
            ylabel('NT Vesicles (1)')

            legend({ sprintf('Inter-peak-times Cleft (%s)', SR{i}), sprintf('Inter-peak-times Nerve Volt (%s)', SR{i}) });

            if runopt.save_figures
                mySaveFig( hfig(end), ...
                    sprintf('InterTimesHistogram_%s',SR{i}), ...
                    runopt.path.nerve, ...
                    runopt.path.nerve, 'width', '0.6\textwidth' );
            end
        end
    
    end
end


%% Phase-based statistics

% [PHI, I, VOLTAGE, SIGNAL] = loadv(src, ...
%     'PHI', 'I', 'VOLTAGE', 'SIGNAL');
% 
% for i = 1:length(SR)
%     if any(strcmp(SR{i}, runopt.SRsToPlot)) || any(strcmp('all', runopt.SRsToPlot))
% 
%         % CaConc
%         for k = find(strcmp({CaConc.sr}, SR{i}))
%             hfig(end+1) = myPhasePlot(CaConc(i).phase, VOLTAGE, PHI, I, ...
%                 '',...sprintf('Mean_R cleft PERIOD (%s)', SR{i}), ...
%                 '',...sprintf('Mean_R cleft PERIOD (%s)', SR{i}), ...
%                 'Mean_\Phi Ca Conc (M)', ...
%                 'Mean_\Phi IHC voltage (V)', ...
%                 'Std_\Phi Ca Conc (M)', ...
%                 'Std_\Phi IHC voltage (V)', plotopt);
% 
%             if runopt.save_figures
%                 mySaveFig( hfig(end), ...
%                     sprintf('SynMeanCleftPeriodPhase_%s',SR{i}), ...
%                     runopt.path.synapse, ...
%                     runopt.path.synapse, 'width', '0.3\textwidth', 'height', '0.5\textwidth' );
%             end
%         end
% 
%         % Cleft
%         for k = find(strcmp({Cleft.sr}, SR{i}))
%             hfig(end+1) = myPhasePlot(Cleft(i).phase, VOLTAGE, PHI, I, ...
%                 '',...sprintf('Mean_R cleft PERIOD (%s)', SR{i}), ...
%                 '',...sprintf('Mean_R cleft PERIOD (%s)', SR{i}), ...
%                 'Mean_\Phi NT Vesicles (1)', ...
%                 'Mean_\Phi IHC voltage (V)', ...
%                 'Std_\Phi NT Vesicles (1)', ...
%                 'Std_\Phi IHC voltage (V)', plotopt);
% 
%             if runopt.save_figures
%                 mySaveFig( hfig(end), ...
%                     sprintf('SynMeanCleftPeriodPhase_%s',SR{i}), ...
%                     runopt.path.synapse, ...
%                     runopt.path.synapse, 'width', '0.3\textwidth', 'height', '0.5\textwidth' );
%             end
%         end
% 
%         % Volt
%         for k = find(strcmp({Volt.sr}, SR{i}))
%             hfig(end+1) = myPhasePlot(Volt(i).phase, VOLTAGE, PHI, I, ...
%                 '',...sprintf('Mean_R Nerve PERIOD (%s)', SR{i}), ...
%                 '',...sprintf('Mean_R Nerve PERIOD (%s)', SR{i}), ...
%                 'Mean_\Phi Nerve Volt (mV)', ...
%                 'Mean_\Phi IHC voltage (V)', ...
%                 'Std_\Phi Nerve Volt (mV)', ...
%                 'Std_\Phi IHC voltage (V)', plotopt);
% 
%             if runopt.save_figures
%                 mySaveFig( hfig(end), ...
%                     sprintf('NerveMeanVoltPeriodPhase_%s',SR{i}), ...
%                     runopt.path.nerve, ...
%                     runopt.path.nerve, 'width', '0.3\textwidth', 'height', '0.5\textwidth' );
%             end
%         end
%     end
% end

%% Time-based statistics

analysisStartTimeIndex = min([length(time),find(runopt.analysisStart <= time,1)]);
warning('Analysis domain striped to < %s , %s >.', time(analysisStartTimeIndex), time(end));

% [TT_Cleft, C_Cleft, J_Cleft, U_Cleft, I_Cleft, ...
%     TT_Volt, C_Volt, J_Volt, U_Volt, I_Volt] = loadv(src, ...
%         'TT_Cleft', 'C_Cleft', 'J_Cleft', 'U_Cleft', 'I_Cleft', ...
%         'TT_Volt', 'C_Volt', 'J_Volt', 'U_Volt', 'I_Volt');

if false
for i = 1:length(SR)
    for j = 1:numel(stimulus)
        if any(strcmp(SR{i}, runopt.SRsToPlot)) || any(strcmp('all', runopt.SRsToPlot))
                        
            k1 = find(strcmp({Cleft.sr}, SR{i}));
            k2 = find(strcmp({Volt.sr}, SR{i}));
                        
            assert(numel(k1) == numel(k2));
        
            for k = 1:numel(k1)
            
                hfig(end+1) = plotopt.figure;
                subplot(2,1,1);

                myTimeIndexPlot( ...
                    Cleft(k1(k)).C{j}, ...
                    Cleft(k1(k)).U{j}, ...
                    Cleft(k1(k)).I{j}, ...
                    Cleft(k1(k)).J{j}, ...
                    Cleft(k1(k)).TT{j}, ...
                    time, ...
                    sprintf('Mean cleft PERIOD (%s) %d', SR{i}, j), ...
                    'Time (ms)', 'Cleft', 'IHC voltage (V)', stimulus(j), ...
                    plotopt);

        %         mySaveFig( hfig(end), ...
        %             sprintf('SynMeanCleftPeriod_%s',SR{i}), ...
        %             runopt.path.synapse, ...
        %             runopt.path.synapse, 'width', '0.7\textwidth' );

                subplot(2,1,2);

                myTimeIndexPlot( ...
                    Volt(k2(k)).C{j}, ...
                    Volt(k2(k)).U{j}, ...
                    Volt(k2(k)).I{j}, ...
                    Volt(k2(k)).J{j}, ...
                    Volt(k2(k)).TT{j}, ...
                    time, ...
                    sprintf('Mean Nerve Volt PERIOD (%s) %d', SR{i}, j), ...
                    'Time (ms)', 'Volt (mV)', 'IHC voltage (V)', stimulus(j), ...
                    plotopt);

        %         mySaveFig( hfig(end), ...
        %             sprintf('NerveMeanVoltPeriod_%s',SR{i}), ...
        %             runopt.path.nerve, ...
        %             runopt.path.nerve, 'width', '0.7\textwidth' );
                if runopt.save_figures
                    mySaveFig( hfig(end), ...
                        sprintf('MeanCleftVoltPeriod_%s_%d',SR{i}, j), ...
                        runopt.path.nerve, ...
                        runopt, 'width', '0.7\textwidth' );
                end
            end
        end
    end
end
end

%% Show figures

arrayfun(@(h) set(h, 'Visible', 'On'), hfig);

end
%% UTILITY FUNCTIONS

    function [ hfig ] = myJoinedFTplot(f, Y1, Y2, LegendText, plotopt )

        [~,xxx] = max(Y1{1});
        
        
        hfig = plotopt.figure;        
        subplot(2,1,1); % Plot single-sided amplitude spectrum.                        
        hold on
        
        for k=1:length(Y1)
            plot(f{k},Y1{k});
        end
        
        xlabel('Frequency (Hz)')
        xlim([0, 4.1*f{1}(xxx)+eps]);
        
        subplot(2,1,2); % Zoom in        
        hold on
        
        Y3 = [Y1,Y2];
        
        for k=1:length(Y3)
            plot(f{k},Y3{k});
        end
        
        legend(LegendText, 'location', 'northeast');

        xlabel('Frequency (Hz)')        
        xlim([0.9*f{1}(xxx),1.1*f{1}(xxx)+2]);
        
    end

    function [ hfig ] = myFTplot(f, YY, fV, VV, LegendText, TitleText, Subsampling, plotopt )

        if nargin >= 7
            f = f(Subsampling);
            YY = YY(Subsampling);
            fV = fV(Subsampling);
            VV = VV(Subsampling);
        end
    
        hfig = plotopt.figure;
        
        [~,xxx] = max(YY);
        
        subplot(2,1,1); % Plot single-sided amplitude spectrum.
        
        title(TitleText);

        s = 1:4.1*xxx; % subsampling and zoom;
        
        plotyy(f(s),YY(s),fV(s),VV(s));
        
%         ax = plotyy(f,YY,fV,VV);                        
%         d = 3*f(xxx)+eps;
%         set(ax, 'xlim', [0 , d]);
        
        xlabel('Frequency (Hz)')

        legend(LegendText, 'location', 'best');

        subplot(2,1,2); % Zoom in
        
        s = (xxx-50):(xxx+50); % subsampling and zoom;
        s = s(s > 0);
        
        plotyy(f(s),YY(s),fV(s),VV(s));
%         ax = plotyy(f,YY,fV,VV);
%         d = 0.04*f(xxx)+eps;    
%         set(ax, 'xlim', [f(xxx) - d , f(xxx) + d]);
        

        xlabel('Frequency (Hz)')        

    end



    function [ hfig ] = myPhasePlot(SIG_1, SIG_2, PHI, I, T1, T2, YL11, YL12, YL21, YL22, plotopt)
        
        hfig = plotopt.figure;    
        subplot(2,1,1); % plot MEAN

        [hAx, ~, ~] = plotyy(PHI, SIG_1.mean(I) , ...
                             PHI, SIG_2.mean(I) );

        title( T1 )
        xlabel('Phase (rad)')        

        ylabel(hAx(1),YL11) % left y-axis
        ylabel(hAx(2),YL12) % right y-axis

        set(hAx, ...
            'XLim', [0,2*pi], ...
            'XTick', [0: pi/2 : 2*pi], ...
            'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'});


        subplot(2,1,2); % plot STD

        [hAx, ~, ~] = plotyy(PHI, SIG_1.std(I) , ...
                             PHI, SIG_2.std(I) );

        title( T2 )
        xlabel('Phase (rad)')        

        ylabel(hAx(1),YL21) % left y-axis
        ylabel(hAx(2),YL22) % right y-axis

        set(hAx, ...
            'XLim', [0,2*pi], ...
            'XTick', [0: pi/2 : 2*pi], ...
            'XTickLabel', {'0', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'});
        
    end

    function myTimeIndexPlot(C, U, I, J, TT, Time, TIT, XL, YL1, YL2, stimulus, plotopt)
    
        if stimulus.frequency == 0

            hAx = plotyy(TT, C ./ J , ...
                         TT, U ./ I );

        else

            hAx = plotyy(Time(1:length(C)).ms, C ./ I , ...
                         Time(1:length(U)).ms, U ./ I );
        end

        title(TIT)
        xlabel(XL)        

        ylabel(hAx(1),YL1) % left y-axis
        ylabel(hAx(2),YL2) % right y-axis
        
    end





