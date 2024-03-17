classdef Channels < handle
    %CHANNELS 
    
    properties
        num
        state
        
        tau
        tau_blocked
        
        topen
        tclose
        tblocked
        mode

        % markov channels -- see stoch_channels_test_3.m
        alpha
        kp0
        S0t
        V0t
        
    end
    
    methods
        function obj = Channels(n, tau, tau_blocked)
            %CHANNELS Construct an instance of this class
            obj.num = n;
            obj.state = repmat('c', n, 1);
            
            obj.tau = tau;
            obj.tau_blocked = tau_blocked;
            
            % obj.topen = zeros(n,1);
            obj.topen = -inf * ones(n,1);
            obj.tclose = zeros(n,1);
            obj.tblocked = zeros(n,1);
            
            obj.mode = ones(n,1);
        end        
    end
end

