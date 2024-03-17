classdef runOpt < Opt
    
    properties
        
        delete_results_to_save_space = false
        do_not_save_results = false
        
        do_not_save_results_mech = false
        
        no_create_only_load     = false
        
        no_create_only_load_oc_mna = false
        no_create_only_load_oc_mna_circuit_ic = false
        no_create_only_load_oc_mna_dae_ic = false

        no_create_only_load_me_mna_dae_ic = false

        no_create_only_load_synapse_ic = false
        
        calculate_only_ant_IC = false

        read_only_cache (1,:) char = ''
        
        purge char
        CodeVersion
%         nerveID
%         runID
        datetime
        t_exec = struct( ...
            'this', [], ...
            'mech',[], ...
            'oc_mna', [], ...
            'synapse', [], ...
            'nerve', [])
        MNAtexec
        texec_nerve
        texec_mech
        
        figureSaveDir char
        SRsToPlot
        analysisStart Time      = Time(0)
        zeroTime Time           = Time(0)
        
        verbose                 = 2
        Debug                   = false
        
        ReconstructResults      = true
        
        draw_gui                = true
        save_figures            = true
        
        waitbarFunctionAvailable
        
        clear_lock_files_if_found (1,1) logical = false
        
        % default_plot_action = feature('ShowFigureWindows')
        
        name = struct( ...
            'stimulus', '')
        
        path = struct( ...
            'mid', '', ...
            'mech', '', ...
            'oc_mna', '', ...
            'oc_mna_analysis', '', ...
            'synapse', '', ...
            'nerve', '')
        
        hash = struct( ...
            'mid', '', ...
            'mech', '', ...
            'oc_mna', '', ...
            'synapse', '', ...
            'nerve', '' )
        
        found = struct( ...
            'mid', false, ...
            'mech', false, ...
            'oc_mna', false, ...
            'synapse', false, ...
            'nerve', false )
        
        do = struct( ...
            'mid', true, ...
            'mech', true, ...
            'oc_mna', true, ...
            'oc_mna_analysis', true, ...
            'synapse', true, ...
            'nerve', true, ...
            'ant_postprocess', true)
        
        recalculate = struct( ...
            'me_mna_dae_ic', false, ...
            'mid', false, ...
            'mid_statistics', false, ...
            'mech', false, ...
            'mech_statistics', false, ...
            'oc_mna_circuit_ic', false, ...
            'oc_mna_dae_ic', false, ...
        	'oc_mna', false, ...
            'oc_mna_statistics', false, ...
            'synapse_ic', false, ...
        	'synapse', false, ...
            'nerve', false, ...
            'replications', false, ...
            'ant_postprocess', false )
        
        plot = struct( ...
            'mid', false, ...
            'mech', false, ...
            'oc_mna', false, ...
            'ant', false, ...
            'synapse', false, ...
            'synapse_avg', false, ...
            'nerve', false, ...
            'ant_postprocess', false, ...
            'DE', false, ...
            'CS', false, ...
            'FSA', false, ...
            'MFA', false, ...
            'JOIN', false )
        
        analysis_results = struct( ...
            'mid', [], ...
            'mech', [], ...
            'oc_mna', [], ...
            'synapse', [], ...
            'nerve', [] )
        
    end
    properties (Hidden)
        definingProperties = {};
    end
    properties (Constant, Hidden)
        requiredParameters = {};
    end
    methods
        % =================================================================
        % SET FUNCTIONS FOR OPTIONAL PARAMETERS
        % -----------------------------------------------------------------
        function set.purge(obj, val)
            Opt.myassert(val, 'char')
            allowed = {'', 'mechanical', 'electrical', 'synapse', 'nerve'};
            assert(ismember(val, allowed), 'Input must be one of the values %s.', Opt.sprint_set(allowed))
            obj.purge = val;
        end
        function set.CodeVersion(obj, val)
            Opt.myassert(val, 'char')
            
            obj.CodeVersion = val;
        end
%         function set.nerveID(obj, val)
%             Opt.myassert(val, 'scalar', 'numeric', 'natural')
%             
%             obj.nerveID = val;
%         end
%         function set.runID(obj, val)
%             Opt.myassert(val, 'scalar', 'numeric', 'natural')
%             
%             obj.runID = val;
%         end
        function set.datetime(obj, val)
            Opt.myassert(val, 'numeric')
            
            obj.datetime = val;
        end
        % -----------------------------------------------------------------
        % TEXEC
        function set.MNAtexec(obj, val)
            Opt.myassert(val, 'numeric')
            
            obj.MNAtexec = val;
        end
        function set.texec_nerve(obj, val)
            Opt.myassert(val, 'numeric')
            
            obj.texec_nerve = val;
        end
        function set.texec_mech(obj, val)
            Opt.myassert(val, 'numeric')
            
            obj.texec_mech = val;
        end
        % -----------------------------------------------------------------
        % PATH, HASH, FOUND
        function set.path(obj, val)
            f = fields(val);
            for i = 1:numel(f)
                Opt.myassert(val.(f{i}), 'char')
            end
                
            obj.path = val;
        end
        function set.hash(obj, val)
            f = fields(val);
            for i = 1:numel(f)
                Opt.myassert(val.(f{i}), 'char')
            end
                
            obj.hash = val;
        end
        function set.found(obj, val)
            f = fields(val);
            for i = 1:numel(f)
                Opt.myassert(val.(f{i}), 'scalar', 'logical_convertible')
                val.(f{i}) = logical(val.(f{i}));
            end
                
            obj.found = val;
        end
        % -----------------------------------------------------------------
        % DO
        function set.do(obj, val)
            f = fields(val);
            for i = 1:numel(f)
                Opt.myassert(val.(f{i}), 'scalar', 'logical_convertible')
                val.(f{i}) = logical(val.(f{i}));
            end
                
            obj.do = val;
        end
        % -----------------------------------------------------------------
        % PLOT
        function set.SRsToPlot(obj, val)
            Opt.myassert(val, 'cell')
            
            obj.SRsToPlot = val;
        end
        % -----------------------------------------------------------------
        % ANALYSIS
        function set.analysisStart(obj, val)
            % Opt.myassert(val, 'scalar', 'numeric')
            Opt.myassert(val, 'scalar', 'Time')
            
            obj.analysisStart = val;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % SET FUNCTIONS FOR DEFAULT PARAMETERS
        % -----------------------------------------------------------------
        function set.verbose(obj, val)
            Opt.myassert(val, 'scalar')
            assert(ismember(val, [0 1 2 3]), 'Input must be one of the values 0, 1, 2 or 3.')
            
            obj.verbose = val;
        end
        function set.Debug(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.Debug = logical(val);
        end
        % -----------------------------------------------------------------
        % RECALCULATE
        function set.recalculate(obj, val)
            f = fields(val);
            for i = 1:numel(f)
                Opt.myassert(val.(f{i}), 'scalar', 'logical_convertible')
                val.(f{i}) = logical(val.(f{i}));
            end
                
            obj.recalculate = val;
        end
        % -----------------------------------------------------------------
        % RECONSTRUCT RESULTS
        function set.ReconstructResults(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.ReconstructResults = logical(val);
        end
        % -----------------------------------------------------------------
        % PLOT        
        function set.plot(obj, val)
            f = fields(val);
            for i = 1:numel(f)
                Opt.myassert(val.(f{i}), 'scalar', 'logical_convertible')
                val.(f{i}) = logical(val.(f{i}));
            end
                
            obj.plot = val;
        end
        % -----------------------------------------------------------------
        % GUI
        function set.draw_gui(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.draw_gui = logical(val);
        end
        % -----------------------------------------------------------------
        % SAVE FIGURES
        function set.save_figures(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.save_figures = logical(val);
        end
        % -----------------------------------------------------------------
        % WAITBAR
        function set.waitbarFunctionAvailable(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.waitbarFunctionAvailable = logical(val);
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % GET FUNCTIONS
        % -----------------------------------------------------------------
        % WAITBAR
        function val = get.waitbarFunctionAvailable(obj)
         
            if isempty(obj.waitbarFunctionAvailable)
                val = ~( usejava('jvm') && ~feature('ShowFigureWindows') );
            else
                val = obj.waitbarFunctionAvailable;
            end
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        
        function obj = runOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
        end
        
    end

end