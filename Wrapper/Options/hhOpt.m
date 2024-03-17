classdef hhOpt < Opt
    properties
        method
        script
        scaleFactor
        
        RepID
        nSR
        slicesToExcite
        
        NumDiv                  = 1
        save_method             = 'matlab_matfile'
        plotNerveSimulation     = false
        VoltageAmplitudeFactor  = 1
        samplingFrequency       = 44100
    end
    properties (Hidden)
        definingProperties = { ...
            'method',  ...
            'scaleFactor',  ...
            'slicesToExcite', ...
            'VoltageAmplitudeFactor',  ...
            'samplingFrequency'  ...
            }
    end
    properties (Constant, Hidden)
        requiredParameters = {'method', 'script', 'scaleFactor'}
    end
    methods
        % =================================================================
        % SET FUNCTIONS FOR PARAMETERS
        % -----------------------------------------------------------------
        function set.method(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'hhfloor', 'hhnint', 'fxnint', 'cw'}), ...
                'Input must be one of the values hhfloor, hhnint, fxnint, cw.')

            obj.method = val;
        end
        function set.script(obj, val)
            Opt.myassert(val, 'char')
            assert(ismember(val, {'ANT'}), ...
                'Input must be one of the values ANT.')

            obj.script = val;
        end
        function set.scaleFactor(obj, val)
            assert( isstruct(val) || (isscalar(val) && isnumeric(val)), ...
                'Value must be either struct, or logical scalar.')
            
            obj.scaleFactor = val;
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
            assert(ismember(val, {'c_posix', 'matlab_fwrite', 'matlab_matfile', 'struct'}), ...
                'Input must be one of the values c_posix, matlab_fwrite, matlab_matfile.')
            
            obj.save_method = val;
        end
        function set.plotNerveSimulation(obj, val)
            assert( islogical(val) || ( ischar(val) || isnumeric(val)), ...
                'Value must be either logical, char, or numeric.')
            
            obj.plotNerveSimulation = val;
        end
        function set.VoltageAmplitudeFactor(obj, val)
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.VoltageAmplitudeFactor = val;
        end
        function set.samplingFrequency(obj, val)
            %Opt.myassert(val, 'scalar', 'numeric', 'positive')
            Opt.myassert(val, 'scalar', 'Frequency', 'positive')
            
            obj.samplingFrequency = val;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        
        function obj = hhOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
        end
    end

end


