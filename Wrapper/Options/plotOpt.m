classdef plotOpt < Opt
    
    properties
        fun % MECH/MNA/LCRB
        
        doplot                      = feature('ShowFigureWindows')
        style                       = 'screen'
        
        figure_params               = {}
        
        subplot                     = true
        
        JoinReplications            = true
        
        default_plt_action          = true
        
        interactive_selection       = false
        
        do = struct( ...
            ...'general', struct( ...
            ...    'stimulus', true ...
            ...), ...
            'middle_ear', struct( ...
                'output', true ...
            ), ...
            'mechanical', struct( ...
                'BMx', true, ...
                'TMx', true, ...
                'BMx_profile', true, ...
                'TMx_profile', true, ...
                'BMv', true, ...
                'TMv', true, ...
                'BMv_profile', true, ...
                'TMv_profile', true, ...
                'TMf', true, ...
                'TMf_profile', true ...
            ), ...
            'oc_electric', struct( ...
                'stimulus', true, ...
                'gui', true, ...
                'MNAsettings', true, ...
                'ApicalCurrentMesh', true, ...
                'BasalCurrentMesh', true, ...
                'CurrentSteadyState', true, ...
                'VoltageMesh', true, ...
                'VoltageSteadyState', true, ...
                'VoltageMesh_IHC', true, ...
                'VoltageMesh_OHC', true, ...
                'VoltageMesh_IHC_rel', true, ...
                'VoltageMesh_OHC_rel', true, ...
                'VoltageSteadyState_HC', true, ...
                'CochlearMicrophonic', true, ...
                'MaximalCrossection', true, ...
                'FourierTransform', true, ...
                'voltage_profile_one_graph_resistors', true, ...
                'voltage_profile_subplot_resistors', true, ...
                'voltage_profile_one_graph_long_resistors', true, ...
                'voltage_profile_subplot_long_resistors', true, ...
                'current_profile_one_graph_long_resistors', true, ...
                'voltage_profile_one_graph_nodes', true, ...
                'voltage_profile_subplot_nodes', true, ...
                'channels_mesh', true, ...
                'channels_mesh_rel', true, ...
                'channels_profile', true, ...
                'channels_profile_rel', true, ...
                'mech_BMx', true, ...
                'mech_TMx', true, ...
                'mech_BMv', true, ...
                'mech_TMv', true, ...
                'mech_profile_BMx', true, ...
                'mech_profile_TMx', true, ...
                'mech_profile_BMv', true, ...
                'mech_profile_TMv', true ...
            ), ...
            'ANT', struct( ...
                'Channels', true, ...
                'Vars1', true, ...
                'Vars2', true, ...
                'Stores', true, ...
                'Cleft', true, ...
                'CleftNerve', true ...
            ), ...
            'nerve', struct( ...
                'nerve_maximal_crossection', true, ...
                'nerve_voltage_mesh', true ...
            ), ...
            'optimization', struct( ...
                'cost', true, ...
                'bestx_plot3', true, ...
                'allx_plot3', true, ...
                'x_hist', true, ...
                'x_density', true, ...
                'cost_surface', true, ...
                'cost_hist', true ...
            ) ...
            )
        
        xAxisValues                 = 'x'
        view                        = [30,10] % [-35,10], ...[20 10], ...
        colormap                    = magma % jet(64), ...
        % centered_colormap           = colMapGen([1,0,0], [0,0,1], 64, 'midCol', [0.975, 0.975, 0.975])
        centered_colormap           = seismic
        ExcitationMesh              = true
        
        
        tspan
        movingAverageWindowLength
        
        subsamplingFactor           = 1
        subsamplingFactorFilt       = 1
        subsamplingFactorStores     = 1
        
        figure_handles
    end
    properties (Hidden)
        definingProperties = {};
    end
    properties (Constant, Hidden)
        requiredParameters = {};
    end
    properties (Dependent)
        surf_options
    end
    methods
        
        function h = figure(obj, varargin)
            h = figure(obj.figure_params{:}, varargin{:});
            
            % store figure handle
            if isempty(obj.figure_handles)
                obj.figure_handles = h;
            else
                obj.figure_handles(end+1) = h;
            end
        end
        
        function hide_figures_while_plotting(obj)
            obj.figure_params = [obj.figure_params, {'Visible', 'off'}];
        end
        
        function set_do_all(obj, val)
            Opt.myassert(val, 'logical_convertible', 'scalar')
            
            s = fieldnames(obj.do);
            for j = 1:numel(s)
                f = fieldnames(obj.do.(s{j}));
                for i = 1:numel(f)
                    obj.do.(s{j}).(f{i}) = logical(val);
                end
            end
        end
        function surf_options = get.surf_options(obj)
            switch obj.style
                case 'print'
                    surf_options = { ...
                        'EdgeAlpha','0.5', ...
                        'MeshStyle','column', ...
                        'LineWidth',0.25, ...
                        'FaceColor','flat' ... 'interp' ...
                        };
                case 'screen'
                    surf_options = { ...
                        'EdgeColor','none', ...
                        'FaceColor','flat' ... 'interp' ...
                        };
                otherwise
                    error('Unrecognised option plotopt.style = %s', obj.style)
            end
        end
        % =================================================================
        % SET FUNCTIONS FOR PARAMETERS
        % -----------------------------------------------------------------
        % FUN
        function set.fun(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'ME', 'MECH', 'MNA', 'LCRB', 'MNA_optimize_circuit'}), ...
                'Input must be one of the values ME, MECH, MNA or LCRB.')

            obj.fun = val;
        end
        % -----------------------------------------------------------------
        % DOPLOT
        function set.doplot(obj, val)
            Opt.myassert(val, 'logical_convertible', 'scalar')

            obj.doplot = val;
        end
        % -----------------------------------------------------------------
        % STYLE
        function set.style(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'screen', 'print'}), ...
                'Input must be one of the values screen, print.')

            obj.style = val;
        end
        % -----------------------------------------------------------------
        % DEFAULT_PLT_ACTION
        function set.default_plt_action(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.default_plt_action = logical(val);
        end
        % -----------------------------------------------------------------
        % SUBPLOT
        function set.subplot(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.subplot = logical(val);
        end
        % -----------------------------------------------------------------
        % DO
%         function set.do(obj, val)
%             f = fields(val);
%             for i = 1:numel(f)
%                 Opt.myassert(val.(f{i}), 'scalar', 'logical_convertible')
%                 val.(f{i}) = logical(val.(f{i}));
%             end
%                 
%             obj.do = val;
%         end
        % -----------------------------------------------------------------
        % X AXIS VALUES
        function set.xAxisValues(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'x', 'stacks'}), ...
                'Input must be one of the values x or stacks.')

            obj.xAxisValues = val;
        end
        % -----------------------------------------------------------------
        % VIEW
        function set.view(obj, val)
            Opt.myassert(val, 'numeric')
            
            obj.view = val;
        end
        % -----------------------------------------------------------------
        % COLORMAP
        function set.colormap(obj, val)
            Opt.myassert(val, 'numeric')
            
            obj.colormap = val;
        end
        % -----------------------------------------------------------------
        % MOVING AVERAGE FILTER
        function set.movingAverageWindowLength(obj, val)
            Opt.myassert(val, 'Time')
            
            obj.movingAverageWindowLength = val;
        end
        % -----------------------------------------------------------------
        % TSPAN
        function set.tspan(obj, val)
            if ~isempty(val)
                %Opt.myassert(val, 'numeric')
                Opt.myassert(val, 'Time')
            end
            
            obj.tspan = val;
        end
        % -----------------------------------------------------------------
        % SUBSAMPLING
        function set.subsamplingFactor(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.subsamplingFactor = val;
        end
        function set.subsamplingFactorFilt(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.subsamplingFactorFilt = val;
        end
        function set.subsamplingFactorStores(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.subsamplingFactorStores = val;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        function obj = plotOpt(FUN, varargin)
                    
            % gather name-value pair input arguments
            opt = struct('fun', FUN, varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
            % update plotopt.do values with default plot action
            % (true/false)
            obj.set_do_all(obj.default_plt_action)
            
        end
        % =================================================================
        %
        function [T, ts, t0, tf] = time_samples(obj, time, topt, fsampldesired, require_equally_spaced_points)
        % TIME_SAMPLES performs undersampling of time with a frequency
        % defined by fsampledesired and crops the returned time according
        % to plotting options.
        
            if nargin < 5
                require_equally_spaced_points = true;
            end
        
            if isempty(obj.tspan)
                t0 = topt.total.t0;
                tf = topt.total.tf;
                t_ind = 1:length(time.s);
            else
                t0 = obj.tspan(1);
                tf = obj.tspan(2);
                t_ind = find(time >= t0 & time <= tf);
            end

            timems = time(t_ind).ms;

            [~,tind0] = min(abs(timems - t0.ms));
            [~,tindf] = min(abs(timems - tf.ms));

            signalLength = Time(timems(tindf) - timems(tind0), 'ms');

            if ~require_equally_spaced_points
                numberOfPoints = 1 + round(fsampldesired * signalLength);
                T = linspace( timems(tind0), timems(tindf), numberOfPoints ); % time points

                ts = zeros(size(T));
                for i = 1:length(T)

                    t = find(timems >= T(i),1);

                    tt = [inf, timems(t), inf];
                    if t > 1, tt(1) = timems(t-1); end
                    if t < length(T), tt(3) = timems(t+1); end

                    [~,r] = min(abs(tt - T(i)));

                    ts(i) = t + r - 2;
                end
            else
                desiredNumberOfPoints = 1 + round(fsampldesired * signalLength);
                equalSpacing = max([1, floor(length(timems)/desiredNumberOfPoints)]);
                % ts = 1 : equalSpacing : length(timems);
                ts = tind0 : equalSpacing : tindf;

                T = timems(ts);
                T = T';

                %
                % Like this, it may happen, that time(tindf) is not present. No way around
                % it, unfortunately.
            end
        end
        function InteractiveSelection(obj, plotfcn, sections)
            arguments
                obj
                plotfcn
            end
            arguments (Repeating)
                sections
            end
            if isempty(sections)
                sections = fields(obj.do);
            end
            done = false;
            while ~done
                fprintf('Available plots:\n');
                k = 0;
                sec_sizes = zeros(numel(sections),1);
                allf = {};
                for j = 1:numel(sections)
                    section = sections{j};
                    
                    f = fields(obj.do.(section));
                    sec_sizes(j) = numel(f);
                    
                    allf = [allf; f];
                    
                    fprintf('\n%s\n', section);
                
                    for i = 1:numel(f)
                        fprintf('%4d: %s\n', k + i, f{i});
                        obj.do.(section).(f{i}) = false;
                    end
                    k = k + i;
                end
                fprintf('\n')
                str = input('Select plot number(s) or exit (0): ', 's');

                
                n = str2num(str);
                if n == 0
                    done = true;
                else
                    for i = 1:numel(n)
                        section = sections{find(n(i) <= cumsum(sec_sizes),1)};
                        obj.do.(section).(allf{n(i)}) = true;
                    end

                    plotfcn(obj);
                end
            end
        end
    end
    methods (Static)
        function V = subsample(mfile, variable, varargin)
            v = varargin;
            w{numel(varargin)} = '';

            % First check if the call was
            %
            %   subsample(mfile, 'myVariable', ':', ':', ':')
            %
            % which is equivalent to just loading the whole variable
            %
            if all(cellfun(@(x) ischar(x) && strcmp(x,':'), varargin))
                V = load(mfile.Properties.Source, variable);
                V = V.(variable);
                return
            end

            % Otherwise ...
            %

            % We want to check if the indices are equally spaced, because it is
            % possible to index directly into a matfile (without loading it
            % into memmory) only for equally sppaced indices.
            %
            IndicesEquallySpaced = true;
            % check if indices are equally spaced
            for ind = 1:numel(varargin)
                w{ind} = ':';

                if isnumeric(varargin{ind}) && length(varargin{ind}) > 1
                    di = diff(varargin{ind});
                    if any(di ~= di(1))
                        v{ind} = ':';
                        w{ind} = varargin{ind};
                        IndicesEquallySpaced = false;
                    end
                end
            end

            if IndicesEquallySpaced
                fprintf('Indices spaced equally, %s.%s\n', inputname(1), variable);
                V = mfile.(variable)(varargin{:});
            else
                fprintf('Indices NOT spaced equally, %s.%s\n', inputname(1), variable);
                if all(cellfun(@(x) ischar(x) && strcmp(x,':'), v))
                    V = load(mfile.Properties.Source, variable);
                    V = V.(variable);
                else
                    V = mfile.(variable)(v{:});
                end

                V = V(w{:});
            end
        end
        function save_figure(hfig, FOLDER, NAME)
            savefig(hfig, fullfile(FOLDER, NAME));
            saveas(hfig, fullfile(FOLDER, [NAME, '.png']));
        end
        function emph_line_width(src, evnt)

            if evnt.Peer.LineWidth == 3
                evnt.Peer.LineWidth = 0.5;
            else
                evnt.Peer.LineWidth = 3;
            end

        end
    end

end
