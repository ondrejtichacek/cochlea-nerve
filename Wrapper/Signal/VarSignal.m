classdef VarSignal < Signal & SignalWithAmplitude
    %VARSIGNAL 
    
    properties
        name char
        originalAudio double
        originalSamplingFrequency Frequency
        zeroTime Time = Time(0)

        audio double
        audiotime Time
    end
    properties (Dependent)
        main_tspan
    end
    methods (Static)
        function par = getrequiredParameters()
            
            par = [Signal.getrequiredParameters(), { ...
                    'name', ...
                    'originalAudio', ...
                    'originalSamplingFrequency'}];
        end
    end
    methods
        function obj = VarSignal(varargin)
            obj = obj@Signal(varargin{:});

            extendDefiningProperties(obj, { ...
                'amplitude', ...
                'originalAudio', ...
                'originalSamplingFrequency'});

            % resample audio signal to targen sampling frequency
            a = resample(obj.originalAudio, obj.samplingFrequency.Hz, obj.originalSamplingFrequency.Hz);

            a = a(:)';

            a0 = zeros(1, obj.zeroDuration*obj.samplingFrequency);
            a3 = zeros(1, obj.fadeDuration*obj.samplingFrequency);

            a = [a0, a, a3];

            % normalize
            obj.audio = a / max(abs(a));

            % compute corresponding time samples
            NuSa = length(obj.audio);
            obj.audiotime = Time(linspace(0,NuSa * 1/obj.samplingFrequency.Hz, NuSa), 's');
            
            if isempty(obj.tspan)
                obj.tspan = [obj.audiotime(1),obj.audiotime(end)];
            end
        end
        % =================================================================
        % INTERFACES
        function signal = eval(obj, time)
            signal = obj.interpolateAcceleration(time);
        end
        function fun = eval_fcn(obj, time_unit)
            [~, fun] = obj.interpolateAcceleration([], time_unit);
        end
        function tspan = get.main_tspan(obj)
            t0 = obj.audiotime(1);
            tf = obj.audiotime(end);
            
            tspan = [t0, tf];
        end
        % =================================================================
        function t = default_analysis_start_time(obj)
            t = obj.zeroTime;
        end
    end
    methods (Access = 'private')
        function [signal, fun] = interpolateAcceleration(obj, t, time_unit)
            
            if nargin < 3
                time_unit = 'ms';
            end
            
            % get handle for envelope function
            envelope_fun = obj.get_envelope_fun(time_unit);

            T = obj.audiotime.(time_unit);
            S = obj.audio;
            F = griddedInterpolant(T(:), S(:));
            fun = @(t) envelope_fun(t) .* F(t);
            
            if ~isempty(t)
                signal = fun(t.(time_unit));
            else
                signal = [];
            end
            
        end
    end
    
end
