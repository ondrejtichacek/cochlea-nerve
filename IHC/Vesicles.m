classdef Vesicles < handle
    %VESICLES
    
    properties
        num
        state
        close_channels
    end
    properties (SetAccess=private)
        release_events = []
        release_events_max = 0
    end
    
    methods
        function obj = Vesicles(n)
            arguments
                n (1,1) double {mustBePositive, mustBeInteger}
            end
            %VESICLES
            
            obj.num = n;
            obj.state = zeros(n,1);
            obj.close_channels = cell(1,n);
        end
        
        function logReleaseEvent(obj, t)
            
            nt = obj.release_events_max + 1;
            
            if numel(obj.release_events) < nt
                obj.release_events(nt:nt+100) = NaN;
            end
            
            obj.release_events(nt) = t;
           
            obj.release_events_max = obj.release_events_max + 1;
        end
    end
end

