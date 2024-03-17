classdef mechOpt < Opt
    
    properties
        
        % model switch
        identifier char {mustBeMember(identifier,{'Vetesnik', 'Verhulst'})} = 'Vetesnik'
        specifier char {mustBeMember(specifier,{'', 'v0', 'v1', 'v2', 'dev'})} = 'v2'
        
        % discretization
        Numstacks {mustBeInteger, mustBePositive, mustBeScalar}
        samplingFrequency Frequency
        
        % solver parameters
        Solver
        NumDiv = 1
        
        % save parameters
        save_method char {mustBeMember(save_method,{'c_posix', 'matlab_fwrite', 'matlab_matfile'})} ...
            = 'matlab_matfile'
        inter_save_interval
        
        approximation char {mustBeMember(approximation,{'linear', 'nonlinear', 'nonlinear-dampened'})} ...
            = 'nonlinear'
        
        % amplifier = 'none';
        % amplifier = 'mechanic';
        amplifier = 'electric';
        
        % integration = 'standalone'
        integration = 'electric'
        
        prestin_voltage_sensitivity = 1;
        
        % nonlinearity parameters
        % see amplifier_params.m
        nonlinParams = amplifier_params(0.8)
                
        % gain of the cochlear amplifier (spatial parameter)
        lambda = 1.157;

        % gain of the cochlear amplifier (global parameter)
        % gain = 0.87;
        % gain = 0;
        gain = 1;
        
        % multiplicative coefficient for shearing damping force
        shearing_coefficient {mustBeNumeric, mustBeScalar} = 95
        % shearing_coefficient {mustBeNumeric, mustBeScalar} = 25
        
        BM_tuning = @(x) 165.4 * (10.^(2.1*(1-x)) - 0.88); % Greenwood
        % BM_tuning = @(x) 165.4 * (10.^(2.1*(1-x)*0.9) - 0.88);
        % BM_tuning = @(x) 991.395 * (10.^(1.283*(1-x)) - 0.98); % Li et al. 2021 10.1038/s41598-021-83225-w
        
        bmdata_params = {};

        TM_tuning_shift = 0.05;
        TM_damping_coef = 21.50;
        
        % frequency map
        TM_tuning = struct( ...
            'base', 2.4531, ...
            'alpha', -6.36, ...
            'Fo', 21100) ...
            ...'base', 2.2022, ...
            ...'alpha', -7.3094, ...
            ...'Fo', 20899)

        % convert model output to nm
        NMkonst_BM (1,1) double
        NMkonst_OHC_cilia (1,1) double
        
        char_frequency_model_best_fit (1,:) char
        % name of the char frequency model that fits best the data of this
        % model
        

        % -----------------------------------------------------------------
        % OUTER EAR
        % outer_ear_params = struct( ...
        %     'identifier', 'none')
        
        outer_ear_params = struct( ...
            'identifier', 'Meddis', ...
            'args', struct( ...
                'externalResonanceFilters', [ 10 1 1000 4000]));
        
        % -----------------------------------------------------------------
        % MIDDLE EAR
        
        % middle_ear_params = struct( ...
        %     'identifier', 'none', ...
        %     ... 'identifier', 'lopezpoveda', ...
        %     ... 'identifier', 'jepsen' ...
        %     'AMo', 1e-6, ...        % loudnes level
        %     'L', -10);              % shift in dB spl
        
        middle_ear_params = struct( ...
            'identifier', 'PBLL', ...
            ... % voltage to displacement of stapes (oval window)
            'vd_OW', 1e-4)  % [length_unit / V]
        
        
        % store arbitrary info here - it affects the simulation hash
        user_settings = struct()

        % store arbitrary info here - it does not affect the simulation hash
        user_data = struct()
    end
    
    properties (Dependent)
        h
        YCutFun
        YCutFunDer
        algorithm
    end
    
    properties (Transient)
        vihc_ss
        vohc_ss
        vohc_0
        OdeFun
        OdeFun_t
        OdeFun_y
        OdeFun_v
        jacobian
        jacobian_nonlin_part
        jacobian_nonlin_part_vohc
        jacobian_lin_part_vohc
        mass
        xgrid
    end
    
    properties (Constant)
    end
    properties (Dependent, Hidden)
        definingProperties          % see get method
    end
    properties (Constant, Hidden)
        requiredParameters = {};
    end
       
    methods
        % =================================================================
        % SET FUNCTIONS FOR OPTIONAL PARAMETERS
        % -----------------------------------------------------------------
        function set.Solver(obj, val)
            assert(isa(val, 'function_handle') || ischar(val), ...
                'Input must be a function handle or a char, not a `%s`.', class(val))
            
            if ischar(val)
                val = str2func(val);
            end
            
            obj.Solver = val;
        end
        function set.NumDiv(obj, val)
            if ischar(val)
                assert(ismember(val, {'auto'}), ...
                'Input must be one of the values auto.')
            else
                Opt.myassert(val, 'scalar', 'numeric', 'natural')
            end
            
            obj.NumDiv = val;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % GET FUNCTIONS
        % -----------------------------------------------------------------
        function val = get.definingProperties(obj)
            
            val = {...
                'identifier', ...
                'specifier', ...
                'Numstacks', ...
                'samplingFrequency', ...
                'user_settings', ...
                };
            
            switch obj.identifier
                case 'Vetesnik'
                    val = [val, {...
                        'outer_ear_params', ...
                        'middle_ear_params', ...
                        'approximation', ...
                        'integration', ...
                        'amplifier', ...
                        'prestin_voltage_sensitivity', ...
                        'nonlinParams', ...
                        'shearing_coefficient', ...
                        'lambda', ...
                        'gain', ...
                        'NMkonst_BM', ...
                        'NMkonst_OHC_cilia', ...
                    }];
                    switch obj.specifier
                        case 'v1'
                            val = [val, {...
                                'TM_tuning', ...
                            }];
                        case 'v2'
                            val = [val, {...
                                'bmdata_params', ...
                                'BM_tuning', ...
                                'TM_tuning_shift', ...
                                'TM_damping_coef', ...
                            }];
                    end
                otherwise
                    error('identifier %s not recognised', obj.identifier)
            end
        end
        function val = get.YCutFun(obj)
            switch obj.approximation
                case 'linear'
                    val = @(Y) zeros(size(Y));
                case 'nonlinear'
                    val = @(Y) nonlin(Y, ...
                        obj.nonlinParams.y1, obj.nonlinParams.y2, ...
                        obj.nonlinParams.c1, obj.nonlinParams.c2, ...
                        obj.nonlinParams.b, obj.nonlinParams.q);
                case 'nonlinear-dampened'
                    val = @(Y) nonlin(Y, ...
                        obj.nonlinParams.y1, obj.nonlinParams.y2, ...
                        obj.nonlinParams.c1, obj.nonlinParams.c2, ...
                        obj.nonlinParams.b, obj.nonlinParams.q) ...
                        .* (1 - boltzmann(Y*obj.NMkonst, 2, 20 , 2));
            end
        end
        function val = get.YCutFunDer(obj)
            switch obj.approximation
                case 'linear'
                    val = @(Y) zeros(size(Y));
                case 'nonlinear'
                    val = @(Y) nonlin_der(Y, ...
                        obj.nonlinParams.y1, obj.nonlinParams.y2, ...
                        obj.nonlinParams.c1, obj.nonlinParams.c2, ...
                        obj.nonlinParams.b, obj.nonlinParams.q);
                case 'nonlinear-dampened'
                    val = @(Y) nonlin(Y, ...
                        obj.nonlinParams.y1, obj.nonlinParams.y2, ...
                        obj.nonlinParams.c1, obj.nonlinParams.c2, ...
                        obj.nonlinParams.b, obj.nonlinParams.q) ...
                        .* (1 - boltzmann(Y*obj.NMkonst, 2, 20 , 2));
            end
        end
        function val = get.algorithm(obj)
            val = func2str(obj.Solver);
        end
        
        function h = get.h(obj)
            h = 1/obj.samplingFrequency;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        function obj = mechOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
        end
    end
end