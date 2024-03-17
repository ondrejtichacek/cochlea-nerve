classdef timeOpt < Opt
    
    properties
        total
        compute
    end
    properties (Dependent, Hidden)
        definingProperties = {'total'};
    end
    properties (Constant, Hidden)
        requiredParameters = {'total'};
    end
    methods
        % =================================================================
        % SET FUNCTIONS FOR OPTIONAL PARAMETERS
        % -----------------------------------------------------------------
        function set.total(obj, val)
            Opt.myassert(val, 'scalar', 'tspanOpt')
                
            obj.total = val;
        end
        function set.compute(obj, val)
            Opt.myassert(val, 'scalar', 'tspanOpt')
                
            obj.compute = val;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % GET FUNCTIONS
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        % =================================================================
        
        function obj = timeOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
        end
    end
    
end