classdef Signal < stimulusOpt
    %SIGNAL
    
    properties
        Duration Time = Time(0)
        zeroDuration Time = Time(0)
        onsetDuration Time = Time(0)
        offsetDuration Time = Time(0)
        fadeDuration Time = Time(0)
        
        envelope char = 'cos'          % default envelope

        tspan
    end

    properties (Hidden)
        
        definingProperties = { ...
            'samplingFrequency', ...
            'contTime', ...
            'name'};

    end
    properties (Dependent)
        t0
        tf
        fs
    end
    methods (Static)
        function par = getrequiredParameters()
            
            par = [stimulusOpt.getrequiredParameters(), { ...
                    'samplingFrequency'}];
        end
    end
    methods
        function obj = Signal(varargin)
            obj = obj@stimulusOpt(varargin{:});

            extendDefiningProperties(obj, { ...
                'envelope', ...
                'Duration', ...
                'zeroDuration', ...
                'onsetDuration', ...
                'offsetDuration', ...
                'fadeDuration'});

        end
        function t0 = get.t0(obj)
            t0 = obj.tspan(1);
        end
        function tf = get.tf(obj)
            tf = obj.tspan(2);
        end
        function fs = get.fs(obj)
            fs = obj.samplingFrequency;
        end
        % =================================================================
        % PLOT
        function plot(obj, args)
            arguments
                obj
                args.plotAudio = true
                args.plotAcceleration = true
                args.plotSpectrogram = true
            end
            % Plot the signal
            
            numSubplots = 0;
            for i = {'plotAudio', ...
                    'plotAcceleration', ...
                    'plotSpectrogram'}
                if args.(cell2mat(i)) == true
                    numSubplots = numSubplots + 1;
                end
            end
            
            
            % Generate time steps for plotting
            
            if numSubplots > 0
                n = 0;
            
                if args.plotAudio
                    n = n + 1;
                    subplot(numSubplots,1,n)
                    
                    if isempty([obj.audio])
                        warning('Can not plot audio -- original source not available!');
                    else
                        audio = obj(1).audio;
                        for i = 2:numel(obj)
                            audio = audio + obj(i).audio;
                        end
                        plot(obj(1).audiotime.ms, audio, 'Color', 'black')
                        ax = gca;
                        ax.YTick = [];
                        ax.XTick = [];
                        axis tight
                    end
%                     xlabel('Time (ms)');
                    ylabel('Amplitude');
                end
                if args.plotAcceleration
                    n = n + 1;
                    subplot(numSubplots,1,n)
                    
                    if ~isfield(obj, 'acceleration') || isempty([obj.acceleration])
                        warning('Can not plot acceleration -- source not available!');
                    else
                        signal = obj(1).acceleration;
                        for i = 2:numel(obj)
                            signal = signal + obj(i).acceleration;
                        end

                        plot(obj(1).time.ms, signal)
                        xlabel('Time (ms)');
                        ylabel('Acceleration');
                    
                        axis tight
                    end
                end
                if args.plotSpectrogram
                    n = n + 1;
                    subplot(numSubplots,1,n)

                    if isempty([obj.audio])
                        signal = obj(1).acceleration;
                        for i = 2:numel(obj)
                            signal = signal + obj(i).acceleration;
                        end
                    else
                        signal = obj(1).audio;
                        for i = 2:numel(obj)
                            signal = signal + obj(i).audio;
                        end
                    end
                    
                    window = round(numel(signal)/4.5/10);

                    spectrogram(signal,window,[],[],obj(1).samplingFrequency.Hz,'yaxis');
%                     spectrogram(signal,kaiser(32768/2,5),[],32*32768,obj(1).samplingFrequency,'yaxis');
%                     image(t,f,s');
                    ylim([0,20])
                    colorbar off
%                     colorbar('southoutside')
                end
            end
            
        end

        function play(obj, answer)
            % Play the signal audio as sound

            if ~isfield(obj, 'originalAudio') || isempty(obj.originalAudio)
                warning('original source not available!');
                fsaudio = obj.samplingFrequency.Hz;
                s = obj.audio;
            else
                fsaudio = obj.originalSamplingFrequency;
                s = obj.originalAudio;
            end
%                 % resample at 96 kHz
%                 fsaudio = 96000;
%
%                 s = interp1(obj.audiotime, obj.audio, obj.t0 : 1/fsaudio : obj.tf);
                
            
            
            if max(abs(s)) > 0.11
                disp('Warning: Volume of the signal may be too high! Normalize?')
                if nargin < 2
                    answer = input('Enter [Y]/n: ','s');
                end

                switch lower(answer)
                    case {'', 'y', 'yes'}
                        s = 0.11 * s / max(abs(s));
                    case {'n', 'no'}
                        % nothing to do here
                    otherwise
                        error('Unrecognised option `%s`.', answer)
                end
            end

            a = audioplayer(s, fsaudio);
            playblocking(a);
            
        end

        function r = plus(obj1,obj2)
            if all(obj1.tspan == obj2.tspan)
                t = obj1.time;
                
                r = Signal(struct( ...
                    'accelerationsignal', obj1.eval(t) + obj2.eval(t), ...
                    'accelerationtime', t ...
                    ));
                
            else
                error('Different values of `tspan` -- Signals can not be added');
            end
        end
        
    end
    methods (Access = 'protected')
        function fun = get_envelope_fun(obj, time_unit)
            switch obj.envelope
                case 'cos'
                    if isempty(time_unit)
                        fun = @obj.envelope_cos;
                    else
                        [~, t0, t3, e1, e2] = obj.envelope_cos([], time_unit);
                        fun = @(t) Signal.envelope_cos_core(t, t0, t3, e1, e2);
                    end
                case 'click'
                    if isempty(time_unit)
                        fun = @obj.envelope_click;
                    else
                        [~, t0, t3, e1, e2] = obj.envelope_click([], time_unit);
                        fun = @(t) Signal.envelope_cos_core(t, t0, t3, e1, e2);
                    end
                otherwise
                    error('Unknown envelope %s', obj.envelope);
            end
        end
        function [f, t0, t3, e1, e2] = envelope_cos(obj, t, time_unit)
            
            if isa(t, 'Time')
                e1 = obj.onsetDuration;
                e2 = obj.offsetDuration;
                t0 = obj.zeroDuration;
                t3 = obj.tf - obj.fadeDuration;
            else
                e1 = obj.onsetDuration.(time_unit);
                e2 = obj.offsetDuration.(time_unit);
                t0 = obj.zeroDuration.(time_unit);
                t3 = obj.tf.(time_unit) - obj.fadeDuration.(time_unit);
            end
            
            if ~isempty(t)
                f = PureTone.envelope_cos_core(t, t0, t3, e1, e2);
            else
                f = [];
            end
            
            % t_tmp is the time where the envelope reaches 99.9%
            %t_tmp = t0 + e*0.9795;
            %pt =  (1-cos(pi*(t_tmp-t0)/e))/ 2;
            %fprintf('The COS envelope staruates at...')
            %display(t_tmp)
            %plot(t,f,t0,0,'*',e+t0,1,'*',t_tmp,pt,'*')
        end
        function [f, t0, t3, e1, e2] = envelope_click(obj, t, time_unit)
            
            if isa(t, 'Time')
                e1 = obj.onsetDuration;
                e2 = obj.offsetDuration;
                t0 = obj.zeroDuration;
                t3 = t0 + e1 + e2 + obj.clickDuration;
            else
                e1 = obj.onsetDuration.(time_unit);
                e2 = obj.offsetDuration.(time_unit);
                t0 = obj.zeroDuration.(time_unit);
                t3 = t0 + e1 + e2 + obj.clickDuration.(time_unit);
            end
            
            if ~isempty(t)
                f = Signal.envelope_cos_core(t, t0, t3, e1, e2);
            else
                f = [];
            end
        end
    end
    methods (Static)
        function f = envelope_cos_core(t, t0, t3, e1, e2)
            
            t1 = t0 + e1;
            t2 = t3 - e2;

            assert(t2 > t1)
            
            if numel(t) == 1
                if t <= t0
                    f = 0;
                    return
                elseif t > t0 && t < t1
                    f = (1-cos(pi*(t-t0)/e1))/ 2;
                    return
                elseif t >= t1 && t <= t2
                    f = 1;
                    return
                elseif t > t2 && t < t3
                    f = (1+cos(pi*(t-t2)/e2))/ 2;
                    return
                elseif t >= t3
                    f = 0;
                    return
                else
                    error("something is wrong")                
                end
            else            
                f = zeros(size(t));

                % f(t <= t0) = 0; % we don't need to do this, already initialized to zero
                f(t >= t1 & t <= t2) = 1;

                i = (t > t0 & t < t1);

                f(i) = (1-cos(pi*(t(i)-t0)/e1))/ 2;

                i = (t > t2 & t < t3);

                f(i) = (1+cos(pi*(t(i)-t2)/e2))/ 2;

            end
            
        end
    end
    
end
