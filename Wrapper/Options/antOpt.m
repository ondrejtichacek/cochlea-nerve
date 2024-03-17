classdef antOpt < Opt
    properties
        ant
        script

        fiber
        
        RepID
        nSR
        slicesToExcite
        positionsToExcite
        
        NumDiv                      = 1
        save_method                 = 'matlab_matfile'
        numberOfRepetitions         = 1
        plotSynapseSimulation       = true
        plotSynapseVarsAreIdentical = false
        VoltageAmplitudeFactor      = 1
        VoltageAmplitudeAddFactor   = 0
        IgnoreNTRateThreshold       = false
        samplingFrequency           = 44100
        slices                      = {{'max'}} % ! needs to be two "cell" {} brackets as struct function "removes one" to extend dimensions !!!
        transductionopt
        
        force_quick_onset (1,1) logical = false
        
        initialConditions = []
        initialConditions_ss = []

    end
    
    properties (Dependent)
        %fiber
        transductionopt_DefiningData
    end
    
%     properties (Access = private)
%         fiber_private
%     end
    properties (Hidden)
        definingProperties = { ...
            'ant',  ...
            'transductionopt_DefiningData',  ...
            ...'slicesToExcite', ...
            'slices', ...
            'VoltageAmplitudeFactor',  ...
            'VoltageAmplitudeAddFactor',  ...
            'IgnoreNTRateThreshold',  ...
            'force_quick_onset', ...
            'samplingFrequency'  ...
            }
    end
    properties (Constant, Hidden)

        requiredParameters = {'ant', 'script', 'fiber'}
            
    end
    methods(Access = protected)
        % Override copyElement method:
        % https://www.mathworks.com/help/matlab/ref/matlab.mixin.copyable-class.html
        function cpObj = copyElement(obj)
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            
            DeepObjets = {'transductionopt'};
            % Make a deep copy of the DeepCp object
            for i = 1:numel(DeepObjets)
                DeepObj = DeepObjets{i};
                cpObj.(DeepObj) = copy(obj.(DeepObj));
            end
        end
    end
    methods
        
        % =================================================================
        % SET FUNCTIONS FOR PARAMETERS
        % -----------------------------------------------------------------
        function set.ant(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'dev', 'dev_stochastic', 'v3', 'v4', 'eguia', 'meddis', 'sumner', 'sumnerStochastic'}), ...
                'Input must be one of the values dev, eguia, meddis, sumner, sumnerStochastic.')

            obj.ant = val;
        end
        function set.script(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'ANT'}), ...
                'Input must be one of the values ANT.')

            obj.script = val;
        end
        function set.RepID(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.RepID = val;
        end
        function set.nSR(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.nSR = val;
        end
        function set.slicesToExcite(obj, val)
            Opt.myassert(val, 'numeric')
            
            obj.slicesToExcite = val;
        end
        function set.NumDiv(obj, val)
            Opt.myassert(val, 'scalar', 'numeric', 'natural')
            
            obj.NumDiv = val;
        end
        function set.save_method(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'c_posix', 'matlab_fwrite', 'matlab_matfile'}), ...
                'Input must be one of the values c_posix, matlab_fwrite, matlab_matfile.')
            
            obj.save_method = val;
        end
        function set.numberOfRepetitions(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.numberOfRepetitions = val;
        end
        function set.plotSynapseSimulation(obj, val)
            assert( islogical(val) || ( ischar(val) || isnumeric(val)), ...
                'Value must be either logical, char, or numeric.')
            
            obj.plotSynapseSimulation = val;
        end
        function set.plotSynapseVarsAreIdentical(obj, val)
            assert( islogical(val) || ( ischar(val) || isnumeric(val)), ...
                'Value must be either logical, char, or numeric.')
            
            obj.plotSynapseVarsAreIdentical = val;
        end
        function set.VoltageAmplitudeFactor(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.VoltageAmplitudeFactor = val;
        end
        function set.VoltageAmplitudeAddFactor(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.VoltageAmplitudeAddFactor = val;
        end
        function set.IgnoreNTRateThreshold(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.IgnoreNTRateThreshold = val;
        end
        function set.samplingFrequency(obj, val)
            %Opt.myassert(val, 'scalar', 'numeric', 'positive')
            Opt.myassert(val, 'scalar', 'Frequency', 'positive')
            
            obj.samplingFrequency = val;
        end
        function set.slices(obj, val)
            Opt.myassert(val, 'cell')
            
            obj.slices = val;
        end
        
%         function set.fiber(obj, SR)
%             
%             Opt.myassert(SR, 'cell')
%             
%             obj.fiber_private = SR;
%             
%             if numel(SR) > 1
%                 error('Multiple SR not set up')
%             end
%             
%             switch obj.ant
%                 case 'v4'
%                     obj.transductionopt = transductionOpt_v4(obj.ant, SR{1});
%                 otherwise
%                     obj.transductionopt = transductionOpt(obj.ant, SR{1});
%             end
%         end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % GET FUNCTIONS FOR PARAMETERS
        % -----------------------------------------------------------------
        %
%         function val = get.fiber(obj)
%             
%             val = obj.fiber_private;
%             
%         end
        
        function val = get.transductionopt_DefiningData(obj)
            
            val = obj.transductionopt.DefiningData;
            
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        function obj = antOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
        end
    end

end
