classdef mnaOpt < Opt
    properties
        
        capacitance_state_dependence = 'strong'
        
        zero_capacitance_test        = false
        
        IHC_basolateral_conductance_dependance char = 'vihc'
        OHC_basolateral_conductance_dependance char = 'vohc'

        IHC_basolateral_popen_dependance char = 'vihc';
        OHC_basolateral_popen_dependance char = 'vohc';
        
        IHC_MET_dependence char      = 'BM'
        OHC_MET_dependence char      = 'cilia'
        
        ext_basolateral_ihc_channels_popen = true
        ext_basolateral_ohc_channels_popen = true
        
        NumDiv                       = 1
        
        save_method (1,:) char             = 'matlab_matfile' 
        evalAtTimePoints (1,1) logical    = true
        samplingFrequency (1,1) Frequency = Frequency(44.1, 'kHz')
        
        solver char                  = 'ode15s'
        fallback_solver char         = 'ode23t'
        
        solveropt cell               = {}
        
        % tol_ref_spl (1,1) double     = 80;
        tol_ref_spl (1,1) double     = 0;
        volt_tol Voltage             = Voltage(1, 'nV') % 10
        curr_tol Current             = Current(0.1, 'pA') % 1
        
        add_solveropt struct         = struct( ...
                                        'use_jacobian', true, ...
                                        'use_jpattern', false)
        
        num_ohc_per_cs (1,1) double  = 4 % human
                                    
        circuit
        
        IHC_basolateral_channels
        OHC_basolateral_channels
        
        simulation_units = struct( ...
            'time', 's', ...
            'voltage', 'V', ...
            'current', 'A', ...
            'conductance', 'S', ...
            'capacitance', 'F', ...
            'charge', 'C');
        
        OHC_V_rest
        IHC_V_rest
        
        % characteristic frequency
        % cf_model (1,:) char = 'Greenwood 1990'
        cf_model (1,:) char = 'Li 2021'
        cf_parameters = [];
        
        Numstacks
    end
    properties (Dependent)
        channels
        xgrid
        NumstacksScaling
        
    end
    properties (Constant)
    end
    properties (Hidden)
        definingProperties = { ...
            'ext_basolateral_ihc_channels_popen', ...
            'capacitance_state_dependence', ...
            'zero_capacitance_test', ...
            'IHC_basolateral_conductance_dependance', ...
            'OHC_basolateral_conductance_dependance', ...
            'IHC_basolateral_popen_dependance', ...
            'OHC_basolateral_popen_dependance', ...
            'IHC_MET_dependence', ...
            'OHC_MET_dependence', ...
            'Numstacks', ...
            'num_ohc_per_cs', ...
            'circuit', ...
            'IHC_basolateral_channels', ...
            'OHC_basolateral_channels', ...
            'samplingFrequency', ...
            'cf_model', ...
            'cf_parameters', ...
            }
    end
    properties (Constant, Hidden)
        requiredParameters = {'Numstacks'};
    end
    methods        
        % =================================================================
        % GET FUNCTIONS FOR PARAMETERS
        % -----------------------------------------------------------------
        %
        function channels = get.channels(obj)
            channels = [];
            if obj.ext_basolateral_ihc_channels_popen == true
                channels = [channels,  obj.IHC_basolateral_channels];
            end
            if obj.ext_basolateral_ohc_channels_popen == true
                channels = [channels,  obj.OHC_basolateral_channels];
            end
            
            for i = 1:numel(channels)
                channels(i).index = i;
            end
        end
        function val = get.NumstacksScaling(obj)
            val = 3000 / obj.Numstacks;
        end
        function val = get.xgrid(obj)
            % creates an equidistant grid on the x axis
            % val = linspace(0, 1, obj.Numstacks+1);
            % val = val(2:end);
            val = mech.gaussgrid(obj.Numstacks);
            % val = mech.gaussgrid(obj.Numstacks+100);
            % val = val(101:end);
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        
        function obj = mnaOpt(varargin)
            
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
        end
    end


end