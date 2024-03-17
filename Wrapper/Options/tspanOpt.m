classdef tspanOpt < Opt
    
    properties
        t0
        tf
    end
    properties (Dependent)
        tspan
        
        % checkpoint frequency and period are linked by a functional dependence 
        %  -- realized as internal hidden properties
        checkpoint_frequency
        checkpoint_period
    end
    properties (Hidden)
        internal_checkpoint_frequency
        internal_checkpoint_period
    end
    properties (Constant)
    end
    properties (Dependent, Hidden)
        definingProperties = {'t0', 'tf'};
    end
    properties (Constant, Hidden)
        requiredParameters = {'t0', 'tf'};
    end
    
    methods
        % =================================================================
        % SET FUNCTIONS FOR OPTIONAL PARAMETERS
        % -----------------------------------------------------------------
        function set.t0(obj, val)
            %Opt.myassert(val, 'scalar', 'numeric')
            Opt.myassert(val, 'scalar', 'Time')
                
            obj.t0 = val;
        end
        function set.tf(obj, val)
            %Opt.myassert(val, 'scalar', 'numeric')
            Opt.myassert(val, 'scalar', 'Time')
                
            obj.tf = val;
        end
        function set.checkpoint_frequency(obj, val)
            Opt.myassert(val, 'scalar', 'numeric', 'positive')
                
            obj.internal_checkpoint_frequency = val;
            obj.internal_checkpoint_period = Time(1/val, 's');
        end
        function set.checkpoint_period(obj, val)
            %Opt.myassert(val, 'scalar', 'numeric', 'positive')
            Opt.myassert(val, 'scalar', 'Time', 'positive')
                
            obj.internal_checkpoint_frequency = 1/val;
            obj.internal_checkpoint_period = val;
        end
        function set.tspan(obj, val)
            assert(length(val) == 2, 'Length of tspan array must be 2, not %d.', length(val))
            
            obj.t0 = val(1);
            obj.tf = val(2);
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        % =================================================================
        % GET FUNCTIONS
        % -----------------------------------------------------------------
        function tspan = get.tspan(obj)
            tspan = [obj.t0, obj.tf];
        end
        function val = get.checkpoint_frequency(obj)
            
            val = obj.internal_checkpoint_frequency;
        end
        function val = get.checkpoint_period(obj)
            
            val = obj.internal_checkpoint_period;
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        function obj = tspanOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
        end
        
        function [time, ind] = intersect_time(obj, time)
            ind = find(time >= obj.t0 & time <= obj.tf);
            time = time(ind);
        end
        function t = generateTimeSamples(obj, fs, extra_points)
            
            if nargin < 3
                extra_points = 0;
            end
            
            interval = [obj.t0, obj.tf + extra_points * 1/fs];
            
            
            numSteps = ceil(diff(interval) * fs);
            numSamples = numSteps + 1;

            t = linspace(interval(1), interval(2), numSamples);
            
        end
    end
end