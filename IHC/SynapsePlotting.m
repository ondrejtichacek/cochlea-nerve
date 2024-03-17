function SynapsePlotting( IHCVoltage, TIME, ANTStatistics, SynResFile, NerveResFile, antopt, hhopt, opt, runopt, plotopt, slice )
%SYNAPSEPLOTTING

toHz = @(x) x.Hz;

matlab2tikz_options = {...
    'extraaxisoptions', 'scaled ticks=false, tick label style={/pgf/number format/fixed}'};

WindowLength = round(plotopt.movingAverageWindowLength * antopt.samplingFrequency);

if plotopt.subsamplingFactor ~= 1
    warning('Plots subsampled by factor %d', plotopt.subsamplingFactor);
end

%%

if find(strcmp(antopt.ant,{'v3', 'v4', 'dev','dev_stochastic','eguia','sumnerStochastic'}),1) % if eguia OR sumnerStochastic
    
    % init figures
    if plotopt.do.Channels
        fig.Channels        = figure;
    end
    
    if plotopt.JoinReplications
        if plotopt.do.Vars1
            fig.Vars1       = figure;
        end
        if plotopt.do.Vars2
            fig.Vars2       = figure;
        end
        if plotopt.do.Stores
            fig.Stores      = figure;
        end
        if plotopt.do.Cleft
            fig.Cleft       = figure;
        end
        if plotopt.do.CleftNerve
            fig.CleftNerve  = figure;
        end
    else
        fig.Vars1       = [];
        fig.Vars2       = [];
        fig.Stores      = [];
        fig.Cleft       = [];
        fig.CleftNerve  = [];
    end
    
end

if plotopt.do.Channels
    
    if isempty(plotopt.tspan)
        ts = 1;
        te = length(TIME);
    else
        ts = find(plotopt.tspan(1) <= TIME, 1);
        if isempty(ts)
            error('Requested plotting time span [%s, %s] not possible. Select a subset of [%s, %s].', ...
                plotopt.tspan(1), plotopt.tspan(2), TIME(1), TIME(end))
        end

        te = find(plotopt.tspan(2) <= TIME, 1);
        if isempty(te)
            te = length(TIME);
            warning('Requested plotting time span [%s, %s] not possible. Cropping to [%s, %s].', ...
                plotopt.tspan(1), plotopt.tspan(2), plotopt.tspan(1), TIME(end))
        end
    end

    subsampling = unique([ts : plotopt.subsamplingFactor : te, te]);
    subs_filt = unique([ts : plotopt.subsamplingFactorFilt : te, te]);


    % single figures
    % Steady state value of open channels fraction
    % ---------------------------------------------------------------------
    figure(fig.Channels)
    
    cmp = get(groot,'DefaultAxesColorOrder');
    % cmp = rwb(3);
    
    subplot(2,1,1)
    hold on
    plot(TIME(subsampling).ms, IHCVoltage(subsampling,slice), 'Color', cmp(1,:));    
    
    Y = movingAverageFilter( WindowLength, IHCVoltage(:,slice) );
    
    plot(TIME(subs_filt).ms, Y(subs_filt), 'Color', 'k', 'LineWidth', 0.75);
    
    Y = mean(IHCVoltage(:,slice));
    plot(TIME([ts, te]).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
    
    % title('IHC_{ic} voltage')
    ylabel('Voltage (V)')
    xlabel('Time (ms)')
    
    v = abs(max(IHCVoltage(:,slice)) - min(IHCVoltage(:,slice)));
    
    ylim([min(IHCVoltage(:,slice))-0.05*v,max(IHCVoltage(:,slice))+0.05*v+eps]);
    
    % ---------------------------------------------------------------------
    subplot(2,1,2);
    hold on

    vvc = Voltage(min(IHCVoltage(:,slice)) : 0.1 : max(IHCVoltage(:,slice)), 'mV');
    vvl = Voltage(min(-80,min(IHCVoltage(:,slice))-15) : 0.1 : min(IHCVoltage(:,slice)), 'mV');
    vvr = Voltage(max(IHCVoltage(:,slice)) : 0.1 : max(20,max(IHCVoltage(:,slice))+15), 'mV');        
    
    xlim([vvl(1).mV, vvr(end).mV])
    
    Syn = SynResFile{1};
    
    
    plot(vvc.mV, Syn.ChannelsOpenSS{1}(vvc), 'Color', cmp(1,:));
    plot(vvl.mV, Syn.ChannelsOpenSS{1}(vvl), 'Color', cmp(2,:));
    plot(vvr.mV, Syn.ChannelsOpenSS{1}(vvr), 'Color', cmp(3,:));

    % title('Steady state value of open channels fraction')
    ylabel('m_{\infty} (1)')
    xlabel('Voltage (mV)')
    
    if runopt.save_figures
        mySaveFig( fig.Channels, ...
            'volt_and_m_inf', ...
            runopt.path.synapse, ...
            runopt.path.synapse, ...
            'width', '0.4\textwidth', ...
            'height', '0.4\textwidth', ...
            matlab2tikz_options{:});
    end
    
end

%%

tropt = copy(antopt.transductionopt);

switch antopt.ant

    case 'eguia'
        %% transduction EGUIA                     
        
        error('not implemented enough, see sumnerStochastic first');
        
%         for i = 1:length(S) % figures for all SRs
%             
%            
%             switch S(i).name
%                 case 'HSR'
%                     lspec = '-.';
%                 case 'MSR'
%                     lspec = '--';
%                 case 'LSR'
%                     lspec = '-';
%             end
%             
%             
%             %% Vars 1
%             if plotopt.JoinReplications
%                 figure(fig.Vars1);
%                 hold on
%             else
%                 fig.Vars1(i) = figure;
%             end
% 
%             subplot(2,2,2);
%             plot(S(i).T, S(i).m_inf, 'Color', cmp(2,:))
%             title('Steady state value of open channels fraction')
%             ylabel('m_{inf} (1)')
%             xlabel('Time (ms)')
% 
%             subplot(2,2,3);
%             plot(S(i).T, S(i).I, 'Color', cmp(3,:));
%             title('Calcium current density')
%             ylabel('I (mA/cm^{2})') 
%             xlabel('Time (ms)')
% 
%             subplot(2,2,4);
%             plot(S(i).T, S(i).k, 'Color', cmp(4,:), 'LineStyle', lspec);
%             title('Transmitter release function')
%             ylabel('k (1)')
%             xlabel('Time (ms)')
%             
%             %% Vars 2
%             figure(fig.Vars2);
%             hold on
%             subplot(2,1,1)            
%             plot(S(i).T, S(i).m, 'Color', cmp(5,:));
%             ylabel('m (1)')
%             xlabel('Time (ms)')
% 
%             subplot(2,1,2)
%             plot(S(i).T, S(i).C, 'Color', cmp(6,:));
%             ylabel('C (mMol)')
%             xlabel('Time (ms)')
%             
%             %% Stores
%             figure(fig.Stores)
%             hold on
%             plot(S(i).T, S(i).q, lspec);
%             plot(S(i).T, S(i).w, lspec);
%             title('Transduction')
%             xlabel('Time (ms)')
%             
%             leg.Stores{2*i-1} = sprintf('q (%s)', S(i).name);
%             leg.Stores{2*i} = sprintf('w (%s)', S(i).name);
%             
%             %% Cleft
%             figure(fig.Cleft)
%             hold on
%             subplot(length(S),1,i)
% %             plot(S(i).T, S(i).c);
% %             title(sprintf('Cleft (%s)', S(i).name))
% %             xlabel('Time (ms)')
% 
%             VV = S(i).V;
% 
%             [hAx,hLine1,hLine2] = plotyy(S(i).T*1e-3, S(i).c, ...
%                                          S(i).T*1e-3, VV );
% 
%             title(sprintf('Cleft (%s)', S(i).name))
%             xlabel('Time (ms)')
% 
% %             [VV(1), max(VV)*1.1]
%             ylim(hAx(2),[VV(1), max(VV)+abs(max(VV))*0.1])
%             
%             ylabel(hAx(1),'Cleft (c = k*q)') % left y-axis
%             ylabel(hAx(2),'IHC voltage (V)') % right y-axis
%             
%         end
%         
%         if plotopt.JoinReplications
%             for i = 1:length(S)
%                 figure(fig.Stores(i))
%                 legend(leg.Stores((2*i-1):(2*i)),'Location','NorthEastOutside');
%             end
%         else
%             figure(fig.Stores)
%             legend(leg.Stores,'Location','NorthEastOutside');
%         end
        
        

    case {'v3', 'v4', 'dev', 'dev_stochastic', 'sumnerStochastic'}
        
        cmp = get(groot,'DefaultAxesColorOrder');
        % cmp = rwb(6);
        
        if plotopt.JoinReplications
            cmp_stores = rwb(2 * numel(SynResFile));
        else
            cmp_stores = rwb(2);
        end
        
        % cmpa = get(groot,'DefaultAxesColorOrder');
        cmap = rwb(2);
        
        
        leg.Stores = [];
        
        for i = 1:numel(SynResFile)
            
            Syn = SynResFile{i}; % matfile handle
            Nerve = NerveResFile{i}; % matfile handle
            
            lspec = '-';
            
            TT = Syn.t;
            TTnerve = Nerve.t;

            if isempty(plotopt.tspan)
                ts = 1;
                te = length(TT);
            else
                ts = find(plotopt.tspan(1) <= TT, 1);
                if isempty(ts)
                    error('Requested plotting time span [%s, %s] not possible. Select a subset of [%s, %s].', ...
                        plotopt.tspan(1), plotopt.tspan(2), TT(1), TT(end))
                end
                
                te = find(plotopt.tspan(2) <= TT, 1);
                if isempty(te)
                    te = length(TT);
                    warning('Requested plotting time span [%s, %s] not possible. Cropping to [%s, %s].', ...
                        plotopt.tspan(1), plotopt.tspan(2), plotopt.tspan(1), TT(end))
                end
            end

            % subsampling = unique([ts : plotopt.subsamplingFactor : te, te]) - ts + 1;
            subs_filt = unique([ts : plotopt.subsamplingFactorFilt : te, te]) - ts + 1;
            TT = TT(ts:te);
            
            subs_edge = [1, numel(TT)];
            
            % -------------------------------------------------------------
            % Vars 1
            if plotopt.do.Vars1
                if ~antopt.plotSynapseVarsAreIdentical || i == 1

                    if plotopt.JoinReplications
                        figure(fig.Vars1);
                        hold on
                    else
                        fig.Vars1(i) = figure;
                    end

                    subplot(3,1,1); hold on
                    [S, T] = intelligentSubsampling( Syn.m_inf(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    plot(T.ms, S, 'Color', cmp(2,:))

                    Y = movingAverageFilter( WindowLength, Syn.m_inf(ts:te,slice) );    
                    plot(TT(subs_filt).ms, Y(subs_filt), 'Color', 'k', 'LineWidth', 0.75);

                    Y = mean(Syn.m_inf(ts:te,slice));
                    plot(TT(subs_edge).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
                    
    %                 title('Steady state value of open channels fraction')
                    ylabel('m_{\infty} (1)')
                    ylim([0,1.05*max(S)+eps]);
                    set(gca, 'XTick', []);
    %                 xlabel('Time (ms)')

                    subplot(3,1,2); hold on
                    % [S, T] = intelligentSubsampling( Syn.m(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    [S, T] = deal(Syn.m(ts:te,slice), TT);
                    plot(T.ms, S, 'Color', cmp(5,:), ...
                        'DisplayName', 'open')

                    Y = movingAverageFilter( WindowLength, Syn.m(ts:te,slice) );    
                    plot(TT(subs_filt).ms, Y(subs_filt), 'Color', 'k', 'LineWidth', 0.75, ...
                        'DisplayName', 'open moving average')

                    % Y = mean(Syn.m(ts:te,slice));
                    % plot(TT(subs_edge).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--', ...
                        % 'DisplayName', 'open mean')
                    
                    
                    % [S, T] = intelligentSubsampling( Syn.CaV13_num_inactivated(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    [S, T] = deal(Syn.CaV13_num_inactivated(ts:te,slice), TT);
                    plot(T.ms, S / tropt.num_CaV13, 'Color', [0.3,0.3,0.3], ...
                        'DisplayName', 'inactivated')

                    % [S, T] = intelligentSubsampling( Syn.CaV13_num_normal(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    [S, T] = deal(Syn.CaV13_num_normal(ts:te,slice), TT);
                    plot(T.ms, S / tropt.num_CaV13, 'Color', cmp(6,:), ...
                        'DisplayName', 'normal')

                    % [S, T] = intelligentSubsampling( Syn.CaV13_num_burst(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    [S, T] = deal(Syn.CaV13_num_burst(ts:te,slice), TT);
                    plot(T.ms, S / tropt.num_CaV13, 'Color', 'b', ...
                        'DisplayName', 'burst')
                    
                    % [S, T] = intelligentSubsampling( Syn.Ca_blocked(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    [S, T] = deal(Syn.Ca_blocked(ts:te,slice), TT);
                    plot(T.ms, S, 'Color', 'r', ...
                        'DisplayName', 'blocked')
                    
                    ylabel('m (1)')
                    % ylim([0,1.05*max(S)+eps]);
                    set(gca, 'XTick', []);
    %                 xlabel('Time (ms)')

                    legend()

                    subplot(3,1,3); hold on

                    % I = Syn.I(ts:te,slice);
                    I = Syn.I_total(ts:te,slice);
                    
                    assert(numel(slice) == 1)
                    
                    switch antopt.ant
                        case {'v4'}
                            scale_constant = 1e9;
                            Ylab = 'I (pA)';
                        case {'v3', 'dev', 'dev_stochastic'}
                            scale_constant = 1e3;
                            Ylab = 'I (mA/cm^2)';
                        case 'sumnerStochastic'                
                            scale_constant = 1e12;
                            Ylab = 'I (pA)';
                    end
                    
                    for ii = 1:size(I,2)
                    
                        % [S, T] = intelligentSubsampling( I(:,ii), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                        [S, T] = deal( I(:,ii), TT);
                        plot(T.ms, S*scale_constant, 'Color', [0.8,0.8,0.8]);
                    end
                    
                    if size(I,2) > 1
                        I = mean(I,2);
                        
                        [S, T] = intelligentSubsampling( I, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                        plot(T.ms, S*scale_constant, 'Color', cmp(3,:));
                    end

                    Y = movingAverageFilter( WindowLength, I );
                    
                    plot(TT(subs_filt).ms, Y(subs_filt)*scale_constant, 'Color', 'k', 'LineWidth', 0.75);

                    Y = mean(I);
                    plot(TT(subs_edge).ms, [Y,Y]*scale_constant, 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');

    %                 title('Calcium current')
                    ylabel(Ylab)
                    ylim([1.05*min(S*scale_constant)-eps,0]);
                    xlabel('Time (ms)')                                


                end
            end
            
            % -------------------------------------------------------------
            % Vars 2
            if plotopt.do.Vars2
                if ~antopt.plotSynapseVarsAreIdentical || i == 1

                    if plotopt.JoinReplications
                        figure(fig.Vars2);
                        hold on
                    else
                        fig.Vars2(i) = figure;
                    end

                    subplot(3,1,1); hold on
                    
                    assert(numel(slice) == 1);
                    
                    C = Syn.C(ts:te,slice);
                    
                    for ii = 1:size(C,2) % over vesicles

                        [ SC, T ] = intelligentSubsampling( C(:,ii), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                        plot(T.ms, SC*1e6, 'Color', [0.8,0.8,0.8]);

                        % tmp = ANTStatistics.CaConc;
                        % p = plot(tmp.Location{i}, tmp.Peak{i}*1e6, 'x', 'Color',  cmp(6,:));

                    end
                    
                    if size(C,2) > 1
                        C = mean(C,2);
                        
                        [ SC, T ] = intelligentSubsampling( C, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                        plot(T.ms, SC*1e6, 'Color', cmp(6,:));
                    end
                       
                    
                    Y = movingAverageFilter( WindowLength, C );    
                    plot(TT(subs_filt).ms, Y(subs_filt)*1e6, 'Color', 'k', 'LineWidth', 0.75);

                    Y = mean(C);
                    plot(TT(subs_edge).ms, [Y,Y]*1e6, 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
                    
                    ylabel('C [\muM]')
                    ylim([0,1.05*max(SC*1e6)+eps]);
                    set(gca, 'XTick', []);
    %                 xlabel('Time (ms)')

                
                    subplot(3,1,2); hold on
                    
                    assert(numel(slice) == 1)
                    
                    k = Syn.k(ts:te,slice);%.Hz;
                    
                    for ii = 1:size(k,2) % over vesicles
                        
                        [S, T] = intelligentSubsampling( k(:,ii), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                        plot(T.ms, S, 'Color', [0.8,0.8,0.8], 'LineStyle', lspec);

                    end
                        

    %                 title('Transmitter release function')
                    ylabel('k (Hz)')
                    ylim([0,1.05*max(S)+eps]);
                    xlabel('Time (ms)')
                    
                    if size(k,2) > 1
                        k = mean(k,2);
                        
                        [S, T] = intelligentSubsampling( k, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                        plot(T.ms, S, 'Color', cmp(4,:));
                    end
                    
                    Y = movingAverageFilter( WindowLength, k );    
                    plot(TT(subs_filt).ms, Y(subs_filt), 'Color', 'k', 'LineWidth', 0.75);

                    Y = mean(k);
                    plot(TT(subs_edge).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');

                    % Rate
                    subplot(3,1,3); hold on

    %                 ind = find(TT > 50,1);
                    ind = 1; % ???

                    SC_red = SC(ind:end);

                    SCrange = max(SC_red) - min(SC_red);
                    SCc = linspace(min(SC_red),max(SC_red),100);
                    dCplot = (SCc(2)-SCc(1));
                    SCl = unique([0 : dCplot : SCc(1), SCc(1)]);
                    SCr = SCc(end) : dCplot : max(1*1e-6,min(SC_red) + SCrange*1.2);

                    plot(SCc*1e6, (Syn.TransmitterRelease{1}(SCc)), 'Color', cmp(1,:));
                    plot(SCl*1e6, (Syn.TransmitterRelease{1}(SCl)), 'Color', cmp(2,:));
                    plot(SCr*1e6, (Syn.TransmitterRelease{1}(SCr)), 'Color', cmp(3,:));

                %     title('Steady state value of open channels fraction')
                    ylabel('k_{\infty} (1)')
                    xlabel('C (\muM)')

                end
            end
            
            % -------------------------------------------------------------
            % Vars 3
            if false
                if ~antopt.plotSynapseVarsAreIdentical || i == 1

                    if plotopt.JoinReplications
                        figure(fig.Vars1);
                        hold on
                    else
                        fig.Vars1(i) = figure;
                    end

                    subplot(3,1,1); hold on
                    [S, T] = intelligentSubsampling( Syn.CaV13_num_inactivated(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    plot(T.ms, S, 'Color', cmp(2,:))

                    Y = movingAverageFilter( WindowLength, Syn.CaV13_num_inactivated(ts:te,slice) );    
                    plot(TT(subs_filt).ms, Y(subs_filt), 'Color', 'k', 'LineWidth', 0.75);

                    Y = mean(Syn.CaV13_num_inactivated(ts:te,slice));
                    plot(TT(subs_edge).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
                    
    %                 title('Steady state value of open channels fraction')
                    ylabel('m_{\infty} (1)')
                    ylim([0,1.05*max(S)+eps]);
                    set(gca, 'XTick', []);
    %                 xlabel('Time (ms)')

                    subplot(3,1,2); hold on
                    [S, T] = intelligentSubsampling( Syn.CaV13_num_normal(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    plot(T.ms, S, 'Color', cmp(5,:));

                    Y = movingAverageFilter( WindowLength, Syn.CaV13_num_normal(ts:te,slice) );    
                    plot(TT(subs_filt).ms, Y(subs_filt), 'Color', 'k', 'LineWidth', 0.75);

                    Y = mean(Syn.CaV13_num_normal(ts:te,slice));
                    plot(TT(subs_edge).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
                    
                    ylabel('m (1)')
                    ylim([0,1.05*max(S)+eps]);
                    set(gca, 'XTick', []);
    %                 xlabel('Time (ms)')

                    subplot(3,1,3); hold on
                    [S, T] = intelligentSubsampling( Syn.CaV13_num_burst(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    plot(T.ms, S, 'Color', cmp(5,:));

                    Y = movingAverageFilter( WindowLength, Syn.CaV13_num_burst(ts:te,slice) );    
                    plot(TT(subs_filt).ms, Y(subs_filt), 'Color', 'k', 'LineWidth', 0.75);

                    Y = mean(Syn.CaV13_num_burst(ts:te,slice));
                    plot(TT(subs_edge).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
                    
                    ylabel('m (1)')
                    ylim([0,1.05*max(S)+eps]);
                    set(gca, 'XTick', []);
    %                 xlabel('Time (ms)')


                end
            end
            
            % -------------------------------------------------------------
            % Stores
            if plotopt.do.Stores
                if plotopt.JoinReplications
                    figure(fig.Stores);
                    hold on
                    
                    c1 = cmp_stores(i,:);
                    c2 = cmp_stores(end-i+1,:);
                else
                    fig.Stores(i) = figure;
                    hold on
                    
                    c1 = cmp_stores(1,:);
                    c2 = cmp_stores(2,:);
                end

%                 [S, T] = intelligentSubsampling( Syn.q(':',slice), Syn.t, plotopt.subsamplingFactorStores, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactorStores), 'MinPeakDistance', Time(0.1, 'ms') );
%                 plot(T.ms, S, ...
%                     'Color', c1, ...
%                     'LineStyle', lspec);
                
                q = Syn.q(ts:te,slice);
                    
                for ii = 1:size(q,2) % over vesicles

                    [S, T] = intelligentSubsampling( q(:,ii), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    plot(T.ms, S, 'Color', [0.8,0.8,0.8], 'LineStyle', lspec);

                end
                
                if size(q,2) > 1
                    q = sum(q,2);

                    [S, T] = intelligentSubsampling( q, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                    plot(T.ms, S, 'Color', c1);
                end

                [S, T] = intelligentSubsampling( Syn.w(':',slice), Syn.t, plotopt.subsamplingFactorStores, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactorStores), 'MinPeakDistance', Time(0.1, 'ms') );
                plot(T.ms, S, ...
                    'Color', c2, ...
                    'LineStyle', lspec);
                
    %             title('NT dynamics')
                xlabel('Time (ms)')
                ylabel('NT Vesicles (1)')

                leg.Stores{2*i-1} = sprintf('q (%s)', Syn.name);
                leg.Stores{2*i} = sprintf('w (%s)', Syn.name);
            end
            
            % -------------------------------------------------------------
            % Cleft
            if plotopt.do.Cleft
                if plotopt.JoinReplications
                    figure(fig.Cleft);
                    hold on

                    subplot(length(SynResFile),1,i)

                else
                    fig.Cleft(i) = figure;
                end                        

    %             plot(S(i).T, S(i).c);
    %             title(sprintf('Cleft (%s)', S(i).name))
    %             xlabel('Time (ms)')

                VV = Syn.V(':',slice);

                [ Sc, Tc ] = intelligentSubsampling( Syn.c(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                [ Sv, Tv ] = intelligentSubsampling( Syn.V(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                
                % Thresholding for plotting
                % Sc( Sc > 0 & Sc < 1e-4 ) = 0;
                
                if any(isnan(Sc)) || any(isnan(Sv)) 
                    error('NAN')
                end
                
                [hAx, hp1, hp2] = myCenteredPlot( ...
                    Tc.ms, Sc, ...
                    Tv.ms, Sv.mV, ...
                    0, Sv(1).mV, ...
                    'Time (ms)', '', '', cmap);
                    ...'Time (ms)', 'NT Vesicles (1)', 'IHC Voltage (mV)', cmap);
                
                hold(hAx(1), 'on')
                tmp = ANTStatistics.Cleft;
                if ~isempty(tmp.Location)
                    plot(hAx(1), tmp.Location{i}, tmp.Peak{i}, 'x', 'Color',  cmap(1,:))
                end
                
                set(hAx(1), 'YTick', [0,1,2]);
                set(hAx(2), 'YTick', [-55,-30,100]);
                
                cleanfigure('pruneText',false)
                
            end
            
            if true
                
                figure
                
                subplot(3,1,1)
                plot(TIME(subsampling).ms, IHCVoltage(subsampling,slice), 'Color', cmp(1,:));    
                
                
                subplot(3,1,2)
                
                k = Syn.k(ts:te,slice);%.Hz;
                q = Syn.q(ts:te,slice);
                
                kk = sum(k, 2);
                kq = sum(k .* q, 2);

                hold on
                
                [S, T] = intelligentSubsampling( kk, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                
                plot(T.ms, S, 'Color', cmp(6,:), ...
                    'DisplayName', sprintf('kk = %g', mean(kk) ))

                nw = round(Time(25,'ms') / mean(diff(TT)));

                plot(TT.ms, movingAverageFilter(nw, kk), 'Color', 'k', ...
                    'DisplayName', sprintf('avg kk'))
                
                [S, T] = intelligentSubsampling( kq, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                plot(T.ms, S, 'Color', cmp(3,:), ...
                    'DisplayName', sprintf('kq = %g',mean(kq) ))
                
                plot(TT.ms, movingAverageFilter(nw, kq), 'Color', [.4,.4,.4], ...
                    'DisplayName', sprintf('avg kq'))

                legend
                
                subplot(3,1,3)
                
                [S, T] = intelligentSubsampling( Syn.c_prot(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                
                plot(T.ms, S, 'Color', cmp(5,:))
                
            end
            
            % -------------------------------------------------------------
            % Cleft + Nerve
            if plotopt.do.CleftNerve
                if plotopt.JoinReplications
                    figure(fig.CleftNerve);
                    hold on

                    subplot(length(SynResFile),1,i)

                else
                    fig.CleftNerve(i) = figure;
                end                        

    %             plot(S(i).T, S(i).c);
    %             title(sprintf('Cleft (%s)', S(i).name))
    %             xlabel('Time (ms)')

                VV = Syn.V(':',slice);

                [ Sc, Tc ] = intelligentSubsampling( Syn.c(ts:te,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                [ Sn, Tn ] = intelligentSubsampling( Nerve.V(ts:te,slice), TTnerve(ts:te), plotopt.subsamplingFactor, 'Npeaks', round(length(TTnerve)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );

                [hAx, h1, h2] = myMirrorPlot( ...
                    Tc.ms, Sc, ...
                    Tn.ms, Sn, ...
                    'Time (ms)', '', '', cmap);
                
                % yyaxis left
                hold(hAx(1), 'on')
                % tmp = ANTStatistics.Cleft;
                % release_events = tmp.Location{i};
                % event_magnitude = tmp.Peak{i};

                release_events = Syn.vesicle_release_events * 1e3; % to ms
                event_magnitude = ones(size(release_events));
                
                plot(hAx(1), release_events, event_magnitude, 'x', 'Color',  cmap(1,:), ...
                    'DisplayName', sprintf('ves. release: %d', numel(release_events)))
                
                % set(hAx(1), 'YTick', [0,1,2]);
                
                % yyaxis right
                hold(hAx(2), 'on')
                tmp = ANTStatistics.Volt;
                plot(hAx(2),tmp.Location{i}, tmp.Peak{i}, 'x', 'Color',  cmap(2,:), ...
                    'DisplayName', sprintf('action potential: %d', numel(tmp.Location{i})))
                
                % set(hAx(2), 'YTick', [0,50,100]);
                
                legend('Location', 'Best')
                
            end
        end
        
        if plotopt.do.Stores
            if plotopt.JoinReplications
                figure(fig.Stores)
                legend(leg.Stores,'Location','NorthEastOutside');            
            else
                for i = 1:length(fig.Stores)
                    figure(fig.Stores(i))
                    legend(leg.Stores((2*i-1):(2*i)),'Location','NorthEast');

                end
            end
        end
        
        % SAVE FIGURES
        
%         warning('Saving figures temporarily disabled');
                        
        if runopt.save_figures
        
            for i = 1:length(fig.Vars1)
                mySaveFig( fig.Vars1(i), ...
                        sprintf('synapse_vars_1_%d_%s',i,Syn.name), ...
                        runopt.path.synapse, ...
                        runopt.path.synapse, ...
                        'width', '0.4\textwidth', ...
                        'height', '0.5\textwidth', ...
                        matlab2tikz_options{:});
            end
            for i = 1:length(fig.Vars2)
                mySaveFig( fig.Vars2(i), ...
                        sprintf('synapse_vars_2_%d_%s',i,Syn.name), ...
                        runopt.path.synapse, ...
                        runopt.path.synapse, ...
                        'width', '0.4\textwidth', ...
                        'height', '0.5\textwidth', ...
                        matlab2tikz_options{:});
            end
            for i = 1:length(fig.Stores)
                mySaveFig( fig.Stores(i), ...
                        sprintf('stores_%d_%s',i,Syn.name), ...
                        runopt.path.synapse, ...
                        runopt.path.synapse, ...
                        'width', '0.4\textwidth', ...
                        'height', '0.1\textwidth', ...
                        matlab2tikz_options{:});
                legend('');
            end
            for i = 1:length(fig.Cleft)
                mySaveFig( fig.Cleft(i), ...
                        sprintf('cleft_%d_%s',i,Syn.name), ...
                        runopt.path.synapse, ...
                        runopt.path.synapse, ...
                        'width', '0.4\textwidth', ...
                        'height', '0.1\textwidth', ...
                        'cleanfigure', false, ...
                        matlab2tikz_options{:});
            end
            for i = 1:length(fig.CleftNerve)
                mySaveFig( fig.CleftNerve(i), ...
                        sprintf('cleftNerve_%d_%s',i,Syn.name), ...
                        runopt.path.synapse, ...
                        runopt.path.synapse, ...
                        'width', '0.4\textwidth', ...
                        'height', '0.1\textwidth', ...
                        matlab2tikz_options{:});
            end
            
        end
                  
    case 'meddis'
        %% transduction MEDDIS 1986,88

        
        figure
        subplot(2,1,1);
        plot(TIME.ms,IHCVoltage(:,slice));
        title('IHC_{ic} voltage; MEDDIS 1986,88')
        ylabel('Voltage (V)')
        xlabel('Time (s)')
        
        subplot(2,1,2);
        plot(T,[k,k_tmp]);
        title('Transmitter release function; MEDDIS 1986,88')
        ylabel('k (1)')
        xlabel('Time (s)')
        
        figure
        plotyy(T, [q,w], T, c);
        title('Transduction; MEDDIS 1986,88')
        xlabel('Time (s)')
        legend('q','w','c',...
                      'Location','NorthEastOutside');
  
    case 'sumner'
        %% transduction SUMNER 2002

        figure
        subplot(2,2,1);
        plot(TIME.ms,IHCVoltage(:,slice),'Color', cmp(1,:));
        title('IHC_{ic} voltage')
        ylabel('Voltage (V)')
        xlabel('Time (s)')

        subplot(2,2,2);
        plot(T, m_inf,'Color', cmp(2,:))
        title('Steady state value of open channels fraction')
        ylabel('m_inf (1)')
        xlabel('Time (s)')

        subplot(2,2,3);
        plot(T,I,'Color', cmp(3,:));
        title('Calcium current')
        ylabel('I (A) ?') 
        xlabel('Time (s)')

        subplot(2,2,4);
        plot(T,[k,k_tmp],'Color', cmp(4,:));
        title('Transmitter release function')
        ylabel('k (1)')
        xlabel('Time (s)')

        % figure
        % plotyy(T, [m,q,c,w], T, C );
        % title('Transduction')
        % xlabel('Time (s)')
        % legend('m','q','c','w','C',...
        %               'Location','NorthEastOutside');

        figure
        subplot(2,1,1)
        plot(T, m, 'Color', cmp(5,:))
        ylabel('m (1)')
        xlabel('Time (s)')
        subplot(2,1,2)
        plot(T, C, 'Color', cmp(6,:));
        ylabel('C (?)')
        xlabel('Time (s)')

        figure
        plotyy(T, [q,w], T, c);
        title('Transduction')
        xlabel('Time (s)')
        legend('q','w','c',...
                      'Location','NorthEastOutside');

        % figure
        % plot(T,q+w+c);
        % title('test')
        % ylabel('q+w+c (1)')
        % xlabel('Time (s)')          
    
    otherwise
        
        error('ANT:InvalidArguments', ['Unsupported option for parametr hhopt.ant: ', antopt.ant]);
        
end

%% SAVE FIGURES

if find(strcmp(antopt.ant,{'eguia','sumnerStochastic'}),1) % if eguia OR sumnerStochastic

    % save figures
    savefig(fig.Channels, fullfile(runopt.path.synapse,'Channels'));
    savefig(fig.Vars1, fullfile(runopt.path.synapse, 'Vars1'));
    savefig(fig.Vars2, fullfile(runopt.path.synapse, 'Vars2'));
    savefig(fig.Stores, fullfile(runopt.path.synapse, 'Stores'));
    savefig(fig.Cleft, fullfile(runopt.path.synapse, 'Cleft'));

end
    
end

