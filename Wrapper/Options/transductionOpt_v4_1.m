classdef transductionOpt_v4_1 < Opt
    properties
        SR
        model
        
        d_nernst

        % Relative tolerance for Ca_conc calculation. The lower the
        % rel_tol is, the longer history we need to calculate.
        %
        Ca_conc_rel_tol (1,1) double = 10e-2
        % Ca_conc_rel_tol (1,1) double = 5e-2
        % Ca_conc_rel_tol (1,1) double = 1e-2
        
        % Considering that near the synapse, diffusion is not obstructed by
        % orgaleles, proteins etc...
        % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1435176/
        Ca_diffusion_coefficient = 5.2e-10; % m^2/s ...  diffusion coefficient
        
        % Ca_diffusion_coefficient = 1e-10; % m^2/s ...  diffusion coefficient        
        % 0.7-1.2 e-10      ... https://doi.org/10.1016/S0006-3495(85)83972-2
        % 5.3 e-10;         ... https://doi.org/10.1016/0143-4160(87)90027-3
        
        % Ca_diffusion_coefficient = 1e-11; % m^2/s ...  diffusion coefficient
        % 1-4 e-11?? https://www.cell.com/fulltext/S0006-3495(01)75996-6
        
        C_Ca_background_factor (1,1) double = 1

        G_Ca
        M
        num_release_sites (1,1) double = 14
        num_inactive_release_sites (1,1) double = 0
        
        num_CaV13 (1,1) double
        tau_CaV13 (1,1)
        tau_CaV13_blocked (1,1)
        
        % manufacture rate
        y = Frequency(10, 'Hz'); % Sumner 2002

        % loss rate
        l = Frequency(1290, 'Hz'); % to roughly match the time-constant of H&H
        % l = Frequency(2580, 'Hz'); % Sumner 2002

        % reprocessing rate
        x = Frequency(66.3, 'Hz'); % Sumner 2002
        
        % re-uptake rate
        r = Frequency(3290, 'Hz'); % to roughly match the time-constant of H&H
        % r = Frequency(6580, 'Hz'); % Sumner 2002
       
        % CaV1.3 channels activation
        channels_open_ss_parameters = struct( ...
            'normal', struct( ...
                'apex', {{ ...
                    Voltage(-23.4, 'mV'), ...   % [V] - IHC voltage - half-activation boltzmann function parameter
                    Voltage(8, 'mV'), ... % [V] - IHC voltage - slope boltzmann function parameter
                    0.0, 1, ...
                    }}, ...
                'base', {{ ...
                    Voltage(-25.9, 'mV'), ...   % [V] - IHC voltage - half-activation boltzmann function parameter
                    Voltage(7.4, 'mV'), ...  % [V] - IHC voltage - slope boltzmann function parameter
                    0.0, 1, ...
                    }} ...
            ), ...
            ...
            'burst', struct( ...
                'apex', {{ ...
                    Voltage(0, 'mV'), ...
                    Voltage(1, 'mV'), ... % [V] - IHC voltage - slope boltzmann function parameter
                    0, 1, }}, ...
                'base', {{ ...
                    Voltage(0, 'mV'), ...
                    Voltage(1, 'mV'), ...  % [V] - IHC voltage - slope boltzmann function parameter
                    0, 1, }} ...
            ) ...
            );
        
        % markov channels -- see stoch_channels_test_3.m
        channels_markov_properties = struct( ...
            'alpha', -133, ...
            'kp0', 5e5)

        % % Zampini 2010 Fig. 4
        % % P o(min) = 0.03; P o(max) = 0.27; V 1/2 = −33.9 mV; S = 5.2 mV.
        % single_CaV13_popen_parameters = {
        %     Voltage(-33.9, 'mV'), ...
        %     ...Voltage(-34-20+2, 'mV'), ... adj OT 28.5.2021
        %     ...Voltage(-34-15, 'mV'), ... adj OT 31.5.2021
        %     Voltage(5.2, 'mV') }
        % single_CaV13_popen = @(V, H, S) 0.03 + boltzmann(V, 0.27 - 0.03, H, S)

        transmitter_release_parameters

        use_CaV13_modes = true
        
        ribbon_synapse_properties

        DisplayName
    end
%     properties (Dependent, Hidden)
%         transmitter_half_activation
%         transmitter_slope
%     end
    properties (Hidden)
        definingProperties = { ...
            'SR', 'model', ...
            'd_nernst', ...
            'Ca_conc_rel_tol', ...
            'Ca_diffusion_coefficient', ...
            'C_Ca_background_factor', ...
            'G_Ca', 'M', ...
            'num_release_sites', ...
            'num_inactive_release_sites', ...
            'num_CaV13', ...
            'tau_CaV13', ...
            'channels_markov_properties', ...
            'y', 'l', 'x', 'r', ...
            'channels_open_ss_parameters', ...
            'transmitter_release_parameters', ...
            'use_CaV13_modes', ...
            'ribbon_synapse_properties'
        };
            
        requiredParameters = {};
        
    end
    methods
        
        function params = generate_transmitter_release_parameters(obj, x_pos)
            
            par_base = obj.transmitter_release_parameters.base;
            par_apex = obj.transmitter_release_parameters.apex;
            
            % linear model
            params = cell(numel(par_base), 1);
            for i = 1:numel(par_base)
                params{i} = par_base{i} + (par_apex{i} - par_base{i}) * x_pos;
            end
        end

        function params = generate_channels_open_ss_parameters(obj, mode, x_pos)
            
            par_base = obj.channels_open_ss_parameters.(mode).base;
            par_apex = obj.channels_open_ss_parameters.(mode).apex;
            
            % linear model
            params = cell(numel(par_base), 1);
            for i = 1:numel(par_base)
                params{i} = par_base{i} + (par_apex{i} - par_base{i}) * x_pos;
            end
        end
        
        function rs = visualize(obj, varargin)
            tmp = namedargs2cell(obj.ribbon_synapse_properties);
        
            rs = RibbonSynapse_v4( tmp{:}, ...
                'plotflag', true, ...
                'num_release_sites', obj.num_release_sites, ...
                'num_channels', obj.num_CaV13, ...
                varargin{:} ...
                    );
        end

        function rs = visualize_density(obj, varargin)
            tmp = namedargs2cell(obj.ribbon_synapse_properties);
        
            rs = RibbonSynapse_v4( tmp{:}, ...
                'plotflag', true, ...
                'plot_histograms', false, ...
                'num_release_sites', obj.num_release_sites, ...
                'num_channels', obj.num_CaV13, ...
                varargin{:} ...
                    );

            n = 100;

            XY = cell(1, n);

            for i = 1:n
                rs = RibbonSynapse_v4( tmp{:}, ...
                'plotflag', false, ...
                'num_release_sites', obj.num_release_sites, ...
                'num_channels', obj.num_CaV13, ...
                varargin{:} ...
                    );
                XY{i} = rs.xc(:,1:2);
            end

            XY = cat(1,XY{:});

            X = XY(:,1);
            Y = XY(:,2);

            [N,Xedges,Yedges] = histcounts2(X, Y, 100);
            Xcenters = Xedges(1) + cumsum(diff(Xedges));
            Ycenters = Yedges(1) + cumsum(diff(Yedges));

            colormap(flip(gray(128)))
            imagesc(Xcenters, Ycenters, N', [0, 50*mean(N(:))])
        end

        function plot(obj)
            x_pos = linspace(0,1,6);
            cmap = flipud(rwb(numel(x_pos)));
            
            figure
            hold on
            V = Voltage(linspace(-80, 10), 'mV');
            for i = 1:numel(x_pos)
                params = generate_channels_open_ss_parameters(obj, 'normal', x_pos(i));
                plot(V.mV, CaV13.ChannelsOpenSS(V, params{:}), ...
                    'Color', cmap(i,:), ...
                    'DisplayName', sprintf('x = %.2g', x_pos(i)));
            end
            legend()
            title('Ca_V1.3 open probability -- normal mode')
            subtitle('steady state, at different positions along tonotopy axis')
            xlabel('IHC basolateral transmembrane potential (mV)')
            ylabel('channel openen probability')
            
            
            axes('Position',[.6 .2 .26 .26])
            box on
            hold on
            V = Voltage(linspace(-5, 5), 'mV');
            params = generate_channels_open_ss_parameters(obj, 'burst', 0.5);
            plot(V.mV, CaV13.ChannelsOpenSS(V, params{:}), ...
                'Color', 'k', ...
                'DisplayName', sprintf('x = %.2g', x_pos(i)));
            % for i = 1:numel(x_pos)
            %     params = generate_channels_open_ss_parameters(obj, 'burst', x_pos(i));
            %     plot(V.mV, CaV13.ChannelsOpenSS(V, params{:}), ...
            %         'Color', cmap(i,:), ...
            %         'DisplayName', sprintf('x = %.2g', x_pos(i)));
            % end
            % legend()
            title('burst mode')
            % title('Ca_V1.3 open probability -- burst mode')
            % subtitle('steady state, at different positions along tonotopy axis')
            xlabel('IHC receptor potential (mV)')
            ylabel('channel openen probability')
            
            % -------------------------------------------------------------
            
            figure
            hold on
            V = Voltage(sort(unique([linspace(-10, 50), linspace(-1, 3)])), 'mV');
            
            [p_12, p_21] = CaV13.ModeTransfer(V.V);
            
            plot(V.mV, max(0, p_12), ...
                'Color', 'b', ...
                'DisplayName', 'normal --> burst');
            
            plot(V.mV, p_21, ...
                'Color', 'r', ...
                'DisplayName', 'burst --> normal');
            
            legend()
            title('Ca_V1.3 transfer probability')
            % subtitle('burst --> normal')
            xlabel('IHC receptor potential (mV)')
            ylabel('transfer probability')
            
            %figure
            axes('Position',[.5 .2 .36 .36])
            box on
            hold on
            V = Voltage(sort(unique([linspace(-1, 1, 10), linspace(-0.2, 0.2)])), 'mV');
            
            [p_12, p_21] = CaV13.ModeTransfer(V.V);
            
            plot(V.mV, max(0, p_12), ...
                'Color', 'b', ...
                'DisplayName', 'normal --> burst');
            
            plot(V.mV, p_21, ...
                'Color', 'r', ...
                'DisplayName', 'burst --> normal');
            
            % legend()
            % title('Ca_V1.3 transfer probability')
            % subtitle('burst --> normal')
            xlabel('IHC receptor potential (mV)')
            ylabel('transfer probability')
            
            %

            % figure
            % C = linspace(0, 1000, 256);
            % k = boltzmann(C*1e-6, obj.transmitter_release_parameters{:});

            x_pos = linspace(0,1,6);
            cmap = flipud(rwb(numel(x_pos)));
            
            figure
            hold on
            C = linspace(0, 1000, 256);
            for i = 1:numel(x_pos)
                params = generate_transmitter_release_parameters(obj, x_pos(i));
                plot(C, boltzmann(C*1e-6, params{:}), ...
                    'Color', cmap(i,:), ...
                    'DisplayName', sprintf('x = %.2g', x_pos(i)));
            end
            legend()

            plot(C, k.Hz)
            xlabel('Ca2+ concentration (uM)')
            ylabel('Target release rate (Hz)')

            % -------------------------------------------------------------

            obj.visualize()

        end
        
        function obj = transductionOpt_v4_1(model, fiber_properties)
            arguments
                model = 'v4'
                fiber_properties = struct([]);
            end

            % NOTES: (Zampini 2014)
            %   2 kHz gerbil IHC
            %   NO - Popen steady state = 0.024
            %   G = 15 pS / channel
            %   num CaV13 active zones = 20.6+-1.5
            %   ~ 16000 channels / IHC
            %     => 770 ch /active zone
            %       (180 ch /a z @ basal)
            %   tau_a = 0.77+-0.08 @ -21 mV
            
            obj.G_Ca = 15e-12;             % [S/channel]
            
            obj.tau_CaV13 = Time(500, 'us'); % see Zampini 2010 Table 1
            obj.tau_CaV13_blocked = Time(1, 'ms');
            
            obj.d_nernst = 20e-9;
            % obj.d_nernst = 5e-9;

            obj.transmitter_release_parameters = struct( ...
                'apex', {{ ...
                    Frequency(500, 'Hz'), ...
                    5.0, ...      % Hill's coefficient
                    100e-6}}, ... % Half activation [Ca2+ M]
                'base', {{ ...
                    Frequency(500, 'Hz'), ...
                    5.0, ...
                    100e-6}} ...
                );

            obj.C_Ca_background_factor = 0;

            obj.num_release_sites = 14;         % [#]
            obj.num_CaV13 = 72;        % [#]
            obj.G_Ca = 3 * obj.G_Ca; % FAC !!

            % h1 = [10, 30, 60, 100, 150]/10;
            h1 = 3;

            obj.ribbon_synapse_properties = struct( ...
                'distance_vesicle_membrane', 5, ... nm
                'intervesicle_distance', [60,80], ...
                'channel_radius', 4, ...
                'vesicle_radius', 20 ...
                );

            obj.ribbon_synapse_properties.channel_distribution_method = 'LJ';

            obj.ribbon_synapse_properties.channel_distribution_parameters = ...
                struct('s', 500, ...
                       'epsilon0', 0.1, ...
                       'r0', 130, ...
                       'n0', 4, ...
                       'm1', 10, ...
                       'h1', h1, ...
                       's1', 3);

            obj.transmitter_release_parameters = struct( ...
                'apex', {{ ...
                    Frequency(100, 'Hz'), ...Frequency(500, 'Hz'), ...
                    5.0, ...      % Hill's coefficient
                    80e-6}}, ... % Half activation [Ca2+ M]
                'base', {{ ...
                    Frequency(100, 'Hz'), ...Frequency(500, 'Hz'), ...
                    5.0, ...
                    80e-6}} ...
            );

            %%

            fn = fieldnames(fiber_properties);
            for i = 1:numel(fn)
                obj.(fn{i}) = fiber_properties.(fn{i});
            end
            
            %%

            obj.model = model;
            obj.SR = fiber_properties.SR; % try if set (otherwise error)
            
        end
    end

end