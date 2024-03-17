classdef memOpt < Opt
    properties
        test_mode               = false
        numWorkers              = 1
        globalMemoryLimit       = Inf
        useAllAvailableMemory   = true
        minFreeMem              = 0
        minFreeMemPerWorker     = 0
    end
    properties (Dependent)
        memLimit
        maxFileMem
        maxVarMem
    end
    properties (Hidden, Access = 'private')
        maxFileMem_int          = Inf
        maxVarMem_int           = Inf
    end
    properties (Dependent, Hidden)
        definingProperties = {};
    end
    properties (Constant, Hidden)
        requiredParameters = {};
    end

    methods (Static)
        function mem = parse_from_string(val)
            expression = '(?<value>\d+(?:\.\d+)?)\s*(?<unit>[A-z]+)';
            match = regexp(val, expression, 'names');
            assert(numel(match) == 1)

            % allow also syntax 1M or 1G, meaning 1 MB and 1 GB
            if match.unit(end) ~= 'B'
                match.unit = [match.unit, 'B'];
            end
            
            mem = Quantity(str2double(match.value), match.unit, 'GB', 'B');            
        end
    end
    
    methods
        
        % =================================================================
        % SET FUNCTIONS FOR PARAMETERS
        % -----------------------------------------------------------------
        %
        function set.numWorkers(obj, val)
            if ischar(val)
                val = to_Bytes(val);
            end
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.numWorkers = val;
        end
        function set.globalMemoryLimit(obj, val)
            if ischar(val)
                val = to_Bytes(val);
            end
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.globalMemoryLimit = val;
        end
        function set.useAllAvailableMemory(obj, val)
            Opt.myassert(val, 'scalar', 'logical_convertible')
            
            obj.useAllAvailableMemory = logical(val);
        end
        function set.minFreeMem(obj, val)
            if ischar(val)
                val = to_Bytes(val);
            end
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.minFreeMem = val;
        end
        function set.minFreeMemPerWorker(obj, val)
            if ischar(val)
                val = to_Bytes(val);
            end
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.minFreeMemPerWorker = val;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % SET FUNCTIONS FOR DEPENDENT PARAMETERS
        % -----------------------------------------------------------------
        %        
        function set.maxFileMem(obj, val)
            if ischar(val)
                val = to_Bytes(val);
            end
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.maxFileMem_int = val;
        end
        function set.maxVarMem(obj, val)
            if ischar(val)
                val = to_Bytes(val);
            end
            Opt.myassert(val, 'scalar', 'numeric')
            
            obj.maxVarMem_int = val;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % GET FUNCTIONS FOR DEPENDENT PARAMETERS
        % -----------------------------------------------------------------
        %
        function memLim = get.memLimit(obj)
            
            if obj.test_mode == true
                memLim = 3e9;
                warning('Using test mode, with free memory 3e9 B.')
            else
                memLim = memoryLimit();
            end

            % Lower the free memory if required
            if memLim > obj.minFreeMem
                memLim = memLim - obj.minFreeMem;
            else
                warning('Free memory is less than %d, ignoring memopt.minFreeMem option.', memopt.minFreeMem)
            end
        end
        function m = get.maxFileMem(obj)
            m = min([ ...
                obj.globalMemoryLimit / obj.numWorkers - obj.minFreeMemPerWorker, ...
                obj.maxFileMem_int, ...
                obj.memLimit]);
        end
        function m = get.maxVarMem(obj)
            m = min([ ...
                obj.globalMemoryLimit / obj.numWorkers - obj.minFreeMemPerWorker, ...
                obj.maxFileMem_int, ...
                obj.memLimit]);
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        function obj = memOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
        end
    end
    
end


