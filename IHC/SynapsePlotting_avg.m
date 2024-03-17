function SynapsePlotting_avg( IHCVoltage, TIME, ANTStatistics, SynResFile, NerveResFile, antopt, hhopt, opt, runopt, plotopt, slice )
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

%%

if plotopt.do.Channels
    
    if isempty(plotopt.tspan)
        ts.volt = 1;
        te.volt = length(TIME);
    else
        ts.volt = find(plotopt.tspan(1) <= TIME, 1);
        if isempty(ts.volt)
            error('Requested plotting time span [%s, %s] not possible. Select a subset of [%s, %s].', ...
                plotopt.tspan(1), plotopt.tspan(2), TIME(1), TIME(end))
        end

        te.volt = find(plotopt.tspan(2) <= TIME, 1);
        if isempty(te.volt)
            te.volt = length(TIME);
            warning('Requested plotting time span [%s, %s] not possible. Cropping to [%s, %s].', ...
                plotopt.tspan(1), plotopt.tspan(2), plotopt.tspan(1), TIME(end))
        end
    end

    subsampling = unique([ts.volt : plotopt.subsamplingFactor : te.volt, te.volt]);
    subs_filt.volt = unique([ts.volt : plotopt.subsamplingFactorFilt : te.volt, te.volt]);


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
    
    plot(TIME(subs_filt.volt).ms, Y(subs_filt.volt), 'Color', 'k', 'LineWidth', 0.75);
    
    Y = mean(IHCVoltage(:,slice));
    plot(TIME([ts.volt, te.volt]).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
    
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
        error('not implemented enough, see sumnerStochastic first');                

    case {'v3', 'v4', 'dev', 'dev_stochastic', 'sumnerStochastic'}
        
        cmp = get(groot,'DefaultAxesColorOrder');
        cmp_stores = rwb(2);
        cmap = rwb(2);
        
        leg.Stores = [];
            
        Syn = SynResFile{1}; % matfile handle
        
        lspec = '-';
        
        TT.syn = Syn.t;

        if isempty(NerveResFile)
            parts = {'syn'};
        else
            parts = {'syn', 'nrv'};
            Nerve = NerveResFile{1}; % matfile handle
            TT.nrv = Nerve.t;
        end

        for i = 1:numel(parts)
            part = parts{i};

            if isempty(plotopt.tspan)
                ts.(part) = 1;
                te.(part) = length(TT.(part));
            else
                ts.(part) = find(plotopt.tspan(1) <= TT.(part), 1);
                if isempty(ts.syn)
                    error('Requested plotting time span [%s, %s] not possible. Select a subset of [%s, %s].', ...
                        plotopt.tspan(1), plotopt.tspan(2), TT.(part)(1), TT.(part)(end))
                end
                
                te.(part) = find(plotopt.tspan(2) <= TT.(part), 1);
                if isempty(te.(part))
                    te.(part) = length(TT.(part));
                    warning('Requested plotting time span [%s, %s] not possible. Cropping to [%s, %s].', ...
                        plotopt.tspan(1), plotopt.tspan(2), plotopt.tspan(1), TT.(part)(end))
                end
            end
    
            % subsampling = unique([ts : plotopt.subsamplingFactor : te, te]) - ts + 1;
            subs_filt.(part) = unique([ts.(part) : plotopt.subsamplingFactorFilt : te.(part), te.(part)]) - ts.(part) + 1;
            TT.(part) = TT.(part)(ts.(part):te.(part));
            
            subs_edge.(part) = [1, numel(TT.(part))];
        end

        % -------------------------------------------------------------
        % Vars 1
        if plotopt.do.Vars1                
            fig.Vars1 = figure;

            % .............................................................
            subplot(3,1,1); hold on
           
            S = average_over_replicas(SynResFile, 'm_inf', {ts.syn:te.syn, slice});
            plot(TT.syn.ms, S, 'Color', cmp(2,:))

            Y = movingAverageFilter( WindowLength, S );
            plot(TT.syn(subs_filt.syn).ms, Y(subs_filt.syn), 'Color', 'k', 'LineWidth', 0.75);
            
            Y = mean(S);
            plot(TT.syn(subs_edge.syn).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
            
            ylabel('m_{\infty} (1)')
            ylim([0,1.05*max(S)+eps]);
            set(gca, 'XTick', []);

            % .............................................................
            subplot(3,1,2); hold on

            S = average_over_replicas(SynResFile, 'm', {ts.syn:te.syn, slice});
            plot(TT.syn.ms, S, 'Color', cmp(5,:), ...
                'DisplayName', 'open')

            Y = movingAverageFilter( WindowLength, S );    
            plot(TT.syn(subs_filt.syn).ms, Y(subs_filt.syn), 'Color', 'k', 'LineWidth', 0.75, ...
                'DisplayName', 'open moving average')            
            
            % S = average_over_replicas(SynResFile, 'CaV13_num_inactivated', {ts.syn:te.syn, slice});
            % plot(TT.syn.ms, S / tropt.num_CaV13, 'Color', [0.3,0.3,0.3], ...
            %     'DisplayName', 'inactivated')

            % S = average_over_replicas(SynResFile, 'CaV13_num_normal', {ts.syn:te.syn, slice});
            % plot(TT.syn.ms, S / tropt.num_CaV13, 'Color', cmp(6,:), ...
            %     'DisplayName', 'normal')

            % S = average_over_replicas(SynResFile, 'CaV13_num_burst', {ts.syn:te.syn, slice});
            % plot(TT.syn.ms, S / tropt.num_CaV13, 'Color', 'b', ...
            %     'DisplayName', 'burst')
            
            S = average_over_replicas(SynResFile, 'Ca_blocked', {ts.syn:te.syn, slice});
            plot(TT.syn.ms, S, 'Color', 'r', ...
                'DisplayName', 'blocked')
            
            ylabel('m (1)')
            set(gca, 'XTick', []);

            legend()

            % .............................................................
            subplot(3,1,3); hold on

            % I = Syn.I(ts:te,slice);
            % I = Syn.I_total(ts:te,slice);
            I = average_over_replicas(SynResFile, 'I_total', {ts.syn:te.syn, slice});
            
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
                [S, T] = deal( I(:,ii), TT.syn);
                plot(T.ms, S*scale_constant, 'Color', [0.8,0.8,0.8]);
            end
            
            if size(I,2) > 1
                error('not up to date')
                I = mean(I,2);
                
                % [S, T] = intelligentSubsampling( I, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
                plot(T.ms, S*scale_constant, 'Color', cmp(3,:));
            end

            Y = movingAverageFilter( WindowLength, I );
            
            plot(TT.syn(subs_filt.syn).ms, Y(subs_filt.syn)*scale_constant, 'Color', 'k', 'LineWidth', 0.75);

            Y = mean(I);
            plot(TT.syn(subs_edge.syn).ms, [Y,Y]*scale_constant, 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');

            ylabel(Ylab)
            ylim([1.05*min(S*scale_constant)-eps,0]);
            xlabel('Time (ms)')                                


        end
        
        % -------------------------------------------------------------
        % Vars 2
        if plotopt.do.Vars2

            figure;
            tiledlayout flow
            
            assert(numel(slice) == 1);
            
            [C, Cv] = average_over_replicas(SynResFile, 'C', {ts.syn:te.syn, slice});
            
            for ii = 1:size(C,2) % over vesicles

                nexttile
                % plot(TT.syn.ms, C(:,ii)*1e6);
                shadedErrorBar(TT.syn.ms, C(:,ii)*1e6, Cv(:,ii)*1e6)
            end
            
            ylabel('C [\muM]')
            %ylim([0,1.05*max(SC*1e6)+eps]);
            set(gca, 'XTick', []);

            figure;
            tiledlayout flow
            
            assert(numel(slice) == 1);
            
            [C, Cv] = average_over_replicas(SynResFile, 'q', {ts.syn:te.syn, slice});
            
            for ii = 1:size(C,2) % over vesicles
                nexttile
                % plot(TT.syn.ms, C(:,ii));
                shadedErrorBar(TT.syn.ms, C(:,ii), Cv(:,ii))
            end
            
            ylabel('q []')
            % ylim([0,1.05*max(SC*1e6)+eps]);
            set(gca, 'XTick', []);

            figure;
            tiledlayout flow
            
            assert(numel(slice) == 1);
            
            [C, Cv] = average_over_replicas(SynResFile, 'Cq', {ts.syn:te.syn, slice});
            
            for ii = 1:size(C,2) % over vesicles

                nexttile
                % plot(TT.syn.ms, C(:,ii)*1e6);
                shadedErrorBar(TT.syn.ms, C(:,ii)*1e6, Cv(:,ii)*1e6)
            end
            
            ylabel('C [\muM]')
            % ylim([0,1.05*max(SC*1e6)+eps]);
            set(gca, 'XTick', []);
        
            % .............................................................

            fig.Vars2 = figure;

            % .............................................................
            subplot(3,1,1); hold on
            
            assert(numel(slice) == 1);
            
            C = average_over_replicas(SynResFile, 'C', {ts.syn:te.syn, slice});
            
            for ii = 1:size(C,2) % over vesicles
                SC = C(:,ii);
                plot(TT.syn.ms, SC*1e6, 'Color', [0.8,0.8,0.8]);
            end
            
            if size(C,2) > 1
                C = mean(C,2);
                plot(TT.syn.ms, C*1e6, 'Color', cmp(6,:));
            end
                           
            Y = movingAverageFilter( WindowLength, C );    
            plot(TT.syn(subs_filt.syn).ms, Y(subs_filt.syn)*1e6, 'Color', 'k', 'LineWidth', 0.75);

            Y = mean(C);
            plot(TT.syn(subs_edge.syn).ms, [Y,Y]*1e6, 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');
            
            ylabel('C [\muM]')
            ylim([0,1.05*max(SC*1e6)+eps]);
            set(gca, 'XTick', []);
        
            % .............................................................
            subplot(3,1,2); hold on
            
            assert(numel(slice) == 1)
            
            [k, kv] = average_over_replicas(SynResFile, 'k', {ts.syn:te.syn, slice});
            % k = k.Hz;
            
            for ii = 1:size(k,2) % over vesicles
                % S = k(:,ii);
                % plot(TT.syn.ms, S, 'Color', [0.8,0.8,0.8], 'LineStyle', lspec);
                shadedErrorBar(TT.syn.ms, k(:,ii), kv(:,ii))
            end
                
            ylabel('k (Hz)')
            % ylim([0,1.05*max(S)+eps]);
            xlabel('Time (ms)')
            
            if size(k,2) > 1
                k = mean(k,2);
                plot(TT.syn.ms, k, 'Color', cmp(4,:));
            end
            
            Y = movingAverageFilter( WindowLength, k );    
            plot(TT.syn(subs_filt.syn).ms, Y(subs_filt.syn), 'Color', 'k', 'LineWidth', 0.75);

            Y = mean(k);
            plot(TT.syn(subs_edge.syn).ms, [Y,Y], 'Color', [0.5,0.5,0.5], 'LineWidth', 0.75, 'LineStyle', '--');

            % .............................................................
            subplot(3,1,3); hold on

%                 ind = find(TT > 50,1);
            ind = 1; % ???

            SC = C;

            SC_red = SC(ind:end);

            SCrange = max(SC_red) - min(SC_red);
            SCc = linspace(min(SC_red),max(SC_red),100);
            dCplot = (SCc(2)-SCc(1));
            SCl = unique([0 : dCplot : SCc(1), SCc(1)]);
            SCr = SCc(end) : dCplot : max(1*1e-6,min(SC_red) + SCrange*1.2);

            plot(SCc*1e6, Syn.TransmitterRelease{1}(SCc), 'Color', cmp(1,:));
            plot(SCl*1e6, Syn.TransmitterRelease{1}(SCl), 'Color', cmp(2,:));
            plot(SCr*1e6, Syn.TransmitterRelease{1}(SCr), 'Color', cmp(3,:));

        %     title('Steady state value of open channels fraction')
            ylabel('k_{\infty} (1)')
            xlabel('C (\muM)')

        end
        
        % -------------------------------------------------------------
        % Stores
        if plotopt.do.Stores

            fig.Stores = figure;
            
            c1 = cmp_stores(1,:);
            c2 = cmp_stores(2,:);
            
            xlabel('Time (ms)')

            yyaxis left

            ylabel('free pool (# ves)')

            hold on
            
            % [q, qv] = average_over_replicas(SynResFile, 'q', {ts.syn:te.syn, slice});
            [q, qv] = average_over_replicas_dyn(SynResFile, 'q_dyn', {ts.syn:te.syn, slice});

            for ii = 1:size(q,2) % over vesicles
                plot(TT.syn.ms, q(:,ii), '-', 'Color', [0.8,0.8,0.8], 'DisplayName', '');
            end
            
            if size(q,2) > 1
                q = mean(q,2);
                plot(TT.syn.ms, q, '-', 'Color', c1, 'DisplayName', sprintf('q (%s)', Syn.name));
            end

            yyaxis right

            % [w, wv] = average_over_replicas(SynResFile, 'w', {ts.syn:te.syn, slice});
            [w, wv] = average_over_replicas_dyn(SynResFile, 'w_dyn', {ts.syn:te.syn, slice});

            if size(w,2) > 1
                for ii = 1:size(w,2) % over vesicles
                    plot(TT.syn.ms, w(:,ii), 'Color', [0.2,0.2,0.2], 'DisplayName', '');
                end
            end
            
            w = mean(w,2);
            plot(TT.syn.ms, w, 'Color', c2, 'DisplayName', sprintf('w (%s)', Syn.name));
            
            

            % leg.Stores{2*i-1} = sprintf('q (%s)', Syn.name);
            % leg.Stores{2*i} = sprintf('w (%s)', Syn.name);

            ax = gca();

            ax.YAxis(1).Color = c1;%[1 0 0];
            ax.YAxis(2).Color = c2;%[0 0 1];
        end
        
        % -------------------------------------------------------------
        % Cleft
        if plotopt.do.Cleft
            
            fig.Cleft = figure;

            % Sc = average_over_replicas(SynResFile, 'c', {ts.syn:te.syn, slice});
            Sc = average_over_replicas_dyn(SynResFile, 'c_dyn', {ts.syn:te.syn, slice});
            Sv = average_over_replicas(SynResFile, 'V', {ts.syn:te.syn, slice});
            
            Sv = Voltage(Sv, 'V');

            if any(isnan(Sc)) || any(isnan(Sv)) 
                error('NAN')
            end
            
            [hAx, hp1, hp2] = myCenteredPlot( ...
                TT.syn.ms, Sc, ...
                TT.syn.ms, Sv.mV, ...
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

        if plotopt.do.Cleft
            
            fig.Cleft = figure;
            hold on

            t = Time(4, 'ms');
            ff = 1/t;
            num_replicas = numel(SynResFile);
            n_dyn_rep = SynResFile{1}.n_dyn_rep;

            % ....
            tmp = cellfun(@(Syn) Syn.vesicle_release_events_dyn(), SynResFile, 'UniformOutput', false);
            Sc = [tmp{:}];
            Sc = Sc{1};
            Sc = Sc(:) * 1e3;
                
            [counts, edges] = histcounts( Sc, ...
                'BinWidth', t.ms, ...
                'Normalization', 'count');

            counts = ff.Hz * counts / num_replicas / n_dyn_rep;

            histogram( ...
                'BinEdges', edges, ...
                'BinCounts', counts)%, ...
                ...'FaceColor', 'k', ...
                ...'FaceAlpha', 1)

            % ....
            tmp = ANTStatistics.Volt;
            if ~isempty(tmp)
                tmp = [tmp.Location];
                Sc = cat(1, tmp{:});
                Sc = Sc(:);
                    
                [counts, edges] = histcounts( Sc, ...
                    'BinWidth', t.ms, ...
                    'Normalization', 'count');
    
                counts = ff.Hz * counts / num_replicas / n_dyn_rep;
    
                histogram( ...
                    'BinEdges', edges, ...
                    'BinCounts', counts)%, ...
                    ...'FaceColor', 'k', ...
                    ...'FaceAlpha', 1)
            end
            
            xlabel('Time (ms)');
            ylabel('Frequency (Hz)');

        end
        
        if false %true
            
            figure
            
            subplot(3,1,1)
            plot(TIME(subsampling).ms, IHCVoltage(subsampling,slice), 'Color', cmp(1,:));    
            
            
            subplot(3,1,2)
            
            k = Syn.k(ts.syn:te.syn,slice).Hz;
            q = Syn.q(ts.syn:te.syn,slice);
            
            kk = sum(k, 2);
            kq = sum(k .* q, 2);

            hold on
            
            [S, T] = intelligentSubsampling( kk, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
            
            plot(T.ms, S, 'Color', cmp(6,:), ...
                'DisplayName', sprintf('kk = %g', mean(kk) ))

            nw = round(Time(25,'ms') / mean(diff(TT)));

            plot(TT.syn.ms, movingAverageFilter(nw, kk), 'Color', 'k', ...
                'DisplayName', sprintf('avg kk'))
            
            [S, T] = intelligentSubsampling( kq, TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
            plot(T.ms, S, 'Color', cmp(3,:), ...
                'DisplayName', sprintf('kq = %g',mean(kq) ))
            
            plot(TT.syn.ms, movingAverageFilter(nw, kq), 'Color', [.4,.4,.4], ...
                'DisplayName', sprintf('avg kq'))

            legend
            
            subplot(3,1,3)
            
            [S, T] = intelligentSubsampling( Syn.c_prot(ts.syn:te.syn,slice), TT, plotopt.subsamplingFactor, 'Npeaks', round(length(Syn.t)/plotopt.subsamplingFactor), 'MinPeakDistance', Time(0.1, 'ms') );
            
            plot(T.ms, S, 'Color', cmp(5,:))
            
        end
        
        % -------------------------------------------------------------
        % Cleft + Nerve
        if false % plotopt.do.CleftNerve

            fig.CleftNerve = figure;
            
            Sc = average_over_replicas(SynResFile, 'c', {ts.syn:te.syn, slice});
            Sn = average_over_replicas(NerveResFile, 'V', {ts.nrv:te.nrv, slice});

            [hAx, h1, h2] = myMirrorPlot( ...
                TT.syn.ms, Sc, ...
                TT.nrv.ms, Sn, ...
                'Time (ms)', '', '', cmap);
            
            % yyaxis left
            hold(hAx(1), 'on')
            % tmp = ANTStatistics.Cleft;
            % release_events = tmp.Location{i};
            % event_magnitude = tmp.Peak{i};

            release_events = Syn.vesicle_release_events{1} * 1e3; % to ms
            event_magnitude = ones(size(release_events));
            
            plot(hAx(1), release_events, event_magnitude, 'x', 'Color',  cmap(1,:), ...
                'DisplayName', sprintf('ves. release: %d', numel(release_events)))
            
            % set(hAx(1), 'YTick', [0,1,2]);
            
            % yyaxis right
            hold(hAx(2), 'on')
            % tmp = ANTStatistics.Volt;
            % plot(hAx(2),tmp.Location{i}, tmp.Peak{i}, 'x', 'Color',  cmap(2,:), ...
            %     'DisplayName', sprintf('action potential: %d', numel(tmp.Location{i})))
            
            % set(hAx(2), 'YTick', [0,50,100]);
            
            legend('Location', 'Best')
            
        end
        
        if false % plotopt.do.Stores
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
  
    case 'sumner'
        %% transduction SUMNER 2002   
    
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

function [S, V] = average_over_replicas(SynResFile, method, args)
    tmp = cellfun(@(Syn) Syn.(method)(args{:}), SynResFile, 'UniformOutput', false);

    % tmp fix
    if isa(tmp{1}, 'Frequency')
        tmp = cellfun(@(x) x.Hz, tmp, 'UniformOutput', false);
    elseif isa(tmp{1}, 'Voltage')
        tmp = cellfun(@(x) x.V, tmp, 'UniformOutput', false);
    end

    dim = ndims(tmp{1}) + 1;
    S = cat(dim, tmp{:});
    V = std(S, 0, dim);
    S = mean(S, dim);
end


function [S, V] = average_over_replicas_dyn(SynResFile, method, args)
    tmp = cellfun(@(Syn) Syn.(method)(args{:}), SynResFile, 'UniformOutput', false);
    dim = 1;
    tmp = cat(dim, tmp{:});


    dim = ndims(tmp{1}) + 1;
    S = cat(dim, tmp{:});
    V = std(S, 0, dim);
    S = mean(S, dim);
end
