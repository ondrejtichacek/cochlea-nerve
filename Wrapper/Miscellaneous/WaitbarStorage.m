classdef WaitbarStorage < handle
    %WAITBARSTORAGE 
    
    properties
        t
        old_time
        start_time
        division_index
        waitbar_handle
        
        num_div = 1
        
        update_threshold_time = 60 % [sec]
        update_threshold_progress = 0.01 % [percent/100]
        
        print_threshold_time = 60 % [sec]
        print_threshold_progress = 0.05 %  [percent/100]
    end
    
    methods
        function obj = WaitbarStorage()
        end
    end
end

