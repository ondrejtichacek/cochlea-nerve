classdef solverOpt < Opt
    %SOLVEROPT 
    
    properties
        TimeStep
        OutputFcn % = @(~,~,~) true;
        OutputFcnEvalInterval
        Precompute
        
        tol = struct( ...
            'rel_mean', 1e-6, ... % (relative tolerance)
            'rel_comp', 1e-3 ... % individual components (relative tolerance)
        )
        maxiter = 100
    end
    
    properties (Transient)
        jacobian
        mass
    end
    
    properties (Dependent)
        UseOutputFcn
    end
    properties (Hidden)
        definingProperties = {};    
        
        requiredParameters = { ...
            'TimeStep' ...
            };
    end
    
    methods
        
        % UseOutputFcn
        function val = get.UseOutputFcn(obj)
            if isempty(obj.OutputFcn) || obj.OutputFcnEvalInterval == inf
                val = false;
            else
                val = true;
            end
        end
        
        % OutputFcnEvalInterval
        function val = get.OutputFcnEvalInterval(obj)
            if isempty(obj.OutputFcnEvalInterval)
                val = obj.TimeStep;
            else
                val = obj.OutputFcnEvalInterval;
            end
        end
        
        % Constructor
        function obj = solverOpt(varargin)
            
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
        end
        
    end
    
end

