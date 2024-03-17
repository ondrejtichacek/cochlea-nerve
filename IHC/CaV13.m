classdef CaV13
    %CAV13 
    
    properties
    end
    methods (Static)
        function m_inf = ChannelsOpenSS(V, V_half, S, pmin, pmax)
            %CHANNELSOPENSS 
            %
            %   V       ... input voltage (Receptor potential) [V]
            %   m_inf   ... steady state value of the open channels fraction; m_inf \in <0,1>
            
            m_inf = pmin + boltzmann(V, pmax, V_half, S);
        end
        function [p_12, p_21] = ModeTransfer(Vt, Vt_rel)
            %% V1
            % p_12 = boltzmann_3level_alt(Vt_rel, 1.1, 5.46e-3, 15.46e-3, 2e-3, 10e-3) - 0.1;
            % p_12 = 10*Vt_rel;
            p_12 = 0;
            p_21 = gaussian(Vt_rel, 0.05, 0, 0.05e-3);

            %% V2

%             p_12 = boltzmann(Vt, 1, -25e-3, 10e-3);
%             p_21 = 1 - boltzmann(Vt, 1, -25e-3, 10e-3);
        end
        function p_block = CaProtonBlock(proton)
            p_block = boltzmann(proton, 1, 3, 0.3);
        end


        function compare_popen()

            hfig = figure();
            hold on

            colororder([0 0 1; .5 .5 .5])

            xlabel('Voltage (mV)')
            V = linspace(-80, 0, 256);

            %%
            yyaxis right
            ylabel('Single channel open probability')

            % --- 1 ---
            source = 'Zampini 2010 Fig. 3F';
            H = Voltage(-33.9, 'mV');
            S = Voltage(6.1, 'mV');
            pmin = 0.007;
            pmax = 0.154;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'Color', [.5,.5,.5]);

            % --- 2 ---
            source = 'Zampini 2010 Fig. 4E';
            H = Voltage(-33.9, 'mV');
            S = Voltage(5.2, 'mV');
            pmin = 0.03;
            pmax = 0.27;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'Color', [.5,.5,.5]);

            %%
            yyaxis left
            ylabel('Macroscopic G_{Ca} (norm.)')

            % --- 1 ---
            source = 'Johnson 2008 Fig. 1C (immature, apex)';
            H = Voltage(-34.4, 'mV');
            S = Voltage(8.2, 'mV');
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', '-', ...
                'Marker', 'none', ...
                'Color', linecolors(1));

            % --- 2 ---
            source = 'Johnson 2008 Fig. 1C (immature, base)';
            H = Voltage(-33.1, 'mV');
            S = Voltage(7.6, 'mV');
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', '--', ...
                'Marker', 'none', ...
                'Color', linecolors(1));

            % --- 3 ---
            source = 'Johnson 2008 Fig. 1C (mature, apex)';
            H = Voltage(-23.4, 'mV');
            S = Voltage(8.0, 'mV');
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', ':', ...
                'Marker', 'none', ...
                'Color', linecolors(1));

            % --- 4 ---
            source = 'Johnson 2008 Fig. 1C (mature, base)';
            H = Voltage(-25.9, 'mV');
            S = Voltage(7.4, 'mV');
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', '-.', ...
                'Marker', 'none', ...
                'Color', linecolors(1));

            % --- 5 ---
            source = 'Zampini 2010 Fig. 3F';
            H = Voltage(-34, 'mV');
            S = Voltage(4.5, 'mV');
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', '-', ...
                'Marker', 'none', ...
                'Color', linecolors(2));

            % --- 6 ---
            source = 'Zampini 2014 Fig. 3D';
            H = Voltage(-20.6, 'mV');
            S = Voltage(5.9, 'mV');
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', '-', ...
                'Marker', 'none', ...
                'Color', linecolors(3));

            % --- 7 ---
            source = 'Ohn 2016 Fig. 3A'; % mice
            H = Voltage(-28.1, 'mV'); % +- 1.7
            S = Voltage(7.2, 'mV'); % +- 0.4
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', '-', ...
                'Marker', 'none', ...
                'Color', linecolors(4));

            % --- .. ---
            source = 'this model';
            H = Voltage(-30, 'mV');
            S = Voltage(7, 'mV');
            pmin = 0.0;
            pmax = 1;

            plot(V, CaV13.ChannelsOpenSS(V, H.mV, S.mV, pmin, pmax), ...
                'DisplayName', source, ...
                'LineStyle', '-', ...
                'Marker', 'none', ...
                'LineWidth', 2, ...
                'Color', 'k');

            legend('Location','northwest')

        end
    end
    % ---------------------------------------------------------------------
    methods
        function obj = CaV13()
            %CAV13 
        end
    end
end

