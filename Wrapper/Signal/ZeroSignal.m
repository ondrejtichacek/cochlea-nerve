classdef ZeroSignal < Signal & SignalWithAmplitude
    %ZEROSIGNAL
    
    properties
        % Duration Time = Time(0)
        % zeroDuration Time = Time(0)
        % onsetDuration Time = Time(0)
        % offsetDuration Time = Time(0)
        % fadeDuration Time = Time(0)
    end
    properties (Constant)
        name = 'zero'
    end
    properties (Hidden)
        audio double                   % contains the stimulus data
        audiotime Time                 % contains the stimulus time
    end
    properties (Access = 'private')
        acceleration            % interface of 'stapes acceleration' signal
        accelerationtime        % interface of time samples of 'stapes acceleration' signal
    end
    properties (Dependent)
        signal                  % model input signal = stapes acceleration = second derivative of audio signal
        time                    % time samples of input signal
        main_tspan
    end
    methods (Static)
        function par = getrequiredParameters()
            
            par = [Signal.getrequiredParameters(), { ...
                    'tspan'}];
        end
    end
    methods
        function obj = ZeroSignal(varargin)
            obj = obj@Signal(varargin{:});

            extendDefiningProperties(obj, {});
            
            obj.Duration = obj.tf;

            dt = 1/obj.fs;
            obj.accelerationtime = Time(unique([obj.t0.s : dt.s : obj.tf.s, obj.tf.s]), 's');
            
            obj.acceleration = obj.generateAcceleration(obj.accelerationtime);

            obj.audio = obj.acceleration;
            obj.audiotime = obj.accelerationtime;
        end
        % =================================================================
        % INTERFACES
        function signal = get.signal(obj)
            signal = obj.acceleration;
        end
        function time = get.time(obj)
            time = obj.accelerationtime;
        end
        function tspan = get.main_tspan(obj)
            t0 = obj.zeroDuration + obj.onsetDuration;
            tf = obj.Duration - obj.offsetDuration - obj.fadeDuration;
            
            tspan = [t0, tf];
        end
        function signal = eval(obj, time)
            signal = obj.generateAcceleration(time);
        end
        function fun = eval_fcn(obj, time_unit)
            [~, fun] = obj.generateAcceleration([], time_unit);
        end
        % =================================================================
        function t = default_analysis_start_time(obj)
            t = obj.zeroDuration;
        end
        
    end
    methods (Access = 'private')
        function [signal, fun] = generateAcceleration(obj, t, time_unit)
            
            fun = @(t) zeros(size(t));
            
            signal = fun(t);
            
        end
    end
    
end
