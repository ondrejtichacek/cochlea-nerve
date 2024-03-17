classdef StapesSignal < Signal & SignalWithAmplitude
    %STAPESSIGNAL
    
    properties
        %samplingFrequency Frequency
    end
    properties (Hidden)
        audio                   % interface of 'audio' signal
        audiotime               % interface of time samples of 'audio' signal
        velocity                % interface of 'stapes velocity' signal
        velocitytime            % interface of time samples of 'stapes velocity' signal
        acceleration            % interface of 'stapes acceleration' signal
        accelerationtime        % interface of time samples of 'stapes acceleration' signal
    end
    properties (Dependent)
        signal                  % model input signal = stapes acceleration = second derivative of audio signal
        time                    % time samples of input signal
    end
    methods (Static)
        function par = getrequiredParameters()
            
            par = [Signal.getrequiredParameters(), { ...
                    'accelerationtime', ...
                    'acceleration'}];
        end
    end
    methods
        function obj = StapesSignal(varargin)
            obj = obj@Signal(varargin{:});

            extendDefiningProperties(obj, { ...
                'acceleration'});
            
            if isempty(obj.tspan)
                obj.tspan = [obj.accelerationtime(1),obj.accelerationtime(end)];
            end
        end
        % =================================================================
        % INTERFACES
        function signal = get.signal(obj)
            signal = obj.acceleration;
        end
        function time = get.time(obj)
            time = obj.accelerationtime;
        end
        function signal = eval(obj, time)
            signal = obj.interpolateAcceleration(time);
        end
        function fun = eval_fcn(obj, time_unit)
            [~, fun] = obj.interpolateAcceleration([], time_unit);
        end
        % =================================================================
    end
    methods (Access = 'private')
        function [signal, fun] = interpolateAcceleration(obj, t, time_unit)
            
            if nargin < 3
                time_unit = 'ms';
            end
            
            T = obj.accelerationtime.(time_unit);
            S = obj.acceleration;
            F = griddedInterpolant(T(:), S(:));
            fun = @(t) F(t);
            
            if ~isempty(t)
                signal = fun(t.(time_unit));
            else
                signal = [];
            end
            
        end
    end
    
end
