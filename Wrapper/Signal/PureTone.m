classdef PureTone < Signal & SignalWithAmplitude
    %PURETONE
    
    properties
        waveform char = 'sine'
    end
    properties (Hidden)
        audio double                   % contains the stimulus data
        audiotime Time                 % contains the stimulus time
    end
    properties (Dependent, Hidden, Access = 'private')
        om
    end
    properties (Dependent)
        name
        main_tspan
    end
    methods (Static)
        function par = getrequiredParameters()
            
            par = [Signal.getrequiredParameters(), { ...
                    'tspan', ...
                    'frequency', ...
                    'zeroDuration', ...
                    'onsetDuration', ...
                    'offsetDuration', ...
                    'fadeDuration'}];
        end
    end
    methods
        function obj = PureTone(varargin)
            obj = obj@Signal(varargin{:});

            obj.Duration = obj.tf;

            extendDefiningProperties(obj, { ...
                'waveform', ...
                'amplitude', ...
                'frequency'});

            dt = 1/obj.fs;
            obj.audiotime = Time(unique([obj.t0.s : dt.s : obj.tf.s, obj.tf.s]), 's');
            
            % generate time-series based on audiotime
            obj.audio = obj.generateSignal(obj.audiotime);
        end
        % =================================================================
        % INTERFACES
        function name = get.name(obj)
            name = sprintf('%dHz-%ddBSPL',obj.frequency.Hz, round(obj.amplitude));
        end
        function signal = eval(obj, time)
            signal = obj(1).generateSignal(time);
            for i = 2:numel(obj)
                signal = signal + obj(i).generateSignal(time);
            end
        end
        function fun = eval_fcn(obj, time_unit)
            [~, fun] = obj.generateSignal([], time_unit);
        end
        function fun = envelope_fcn(obj, time_unit)
            [~, ~, fun] = obj.generateSignal([], time_unit);
        end
        function tspan = get.main_tspan(obj)
            t0 = obj.zeroDuration + obj.onsetDuration;
            tf = obj.Duration - obj.offsetDuration - obj.fadeDuration;
            
            tspan = [t0, tf];
        end
        % =================================================================
        function om = get.om(obj)
            om = 2*pi*obj.frequency;
        end
        function t = default_analysis_start_time(obj)
            switch obj.envelope
                case 'cos'
                    t =  obj.zeroDuration(1) + obj.onsetDuration;
                otherwise
                    error('Unknown envelope %s', obj.envelope);
            end
        end
    end
    methods (Access = 'private')
        function [signal, fun, envelope_fun] = generateSignal(obj, t, time_unit)
            if nargin < 3
                time_unit = '';
            end
            
            if isa(t, 'Time')
                omega = obj.om;
            else
                frequency_unit = Time.inverse_unit(time_unit);
                omega = obj.om.(frequency_unit);
            end
            
            % get handle for envelope function
            envelope_fun = obj.get_envelope_fun(time_unit);
            
            % multiply the envelope and the sinusoidal signal
            switch obj.waveform
                case 'sine'
                    fun = @(t) envelope_fun(t) .* sin(omega .* t);
                case 'square'
                    fun = @(t) envelope_fun(t) .* sign(sin(omega .* t));
                % case 'sawtooth'
                otherwise
                    error('unknown value %s', obj.waveform)
            end
            
            signal = fun(t);
        end
    end
end
