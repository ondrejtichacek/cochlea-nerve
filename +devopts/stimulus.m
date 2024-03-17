function [TOPT, signal] = stimulus(GlobalSamplingFrequency, tf_extra, args)
arguments
    GlobalSamplingFrequency (1,1) Frequency
    tf_extra (1,1) % Time or double

    args.strict_signal_length (1,1) logical = false

    args.t0 (1,1) = Time(0)
    args.t2 (1,1) = Time(0)
    
    args.amplitude_unit (1,:) char = 'spl'
    
    args.amplitude = []
    args.frequency = []
    
    args.special (1,:) char = ''
    args.signal_file (1,:) char = ''
    args.signal_var struct = struct()

    % args.onset = 'quick';
    args.onset = 'easy';
    args.offset = 'quick';
    % args.offset = 'easy'

    args.waveform = 'sine'

    args.zeroDuration = Time(0.1, 'ms');
    args.fadeDuration = Time(10, 'ms');

    args.debug (1,1) logical = true

    args.plotflag (1,1) logical = false
end

if isa(args.frequency, 'double')
    args.frequency = Frequency(args.frequency, 'Hz');
end

onset = args.onset;
offset = args.offset;

if ~isempty(args.special)
    onsetDuration = onset;
    offsetDuration = offset;

    % signal_duration = @(tf_extra) tf_extra;
    
else
    onsetDuration = onset_period(onset, args.amplitude, args.frequency);
    offsetDuration = onset_period(offset, args.amplitude, args.frequency);
end

switch args.special
    case {'var', 'file'}
        signal_duration = @(tf_extra) ...
              args.zeroDuration ...
            ...+ onsetDuration ...
            + tf_extra ...
            ...+ offsetDuration ...
            + args.fadeDuration;
    otherwise
        signal_duration = @(tf_extra) ...
              args.zeroDuration ...
            + onsetDuration ...
            + tf_extra ...
            + offsetDuration ...
            + args.fadeDuration;
end


topt_fun = @(tf_extra) timeOpt( ...
    'total', tspanOpt( ...
        't0', args.t0, ...
        'tf', signal_duration(tf_extra), ...
        'checkpoint_period', Time(20, 'ms') ...
        ... 'checkpoint_frequency', Frequency(100, 'Hz') ...
        ) ...
    );

%%

if ~isempty(args.special)
    switch args.special
        case 'zero'
            signal = ...
                ZeroSignal( ...
                    'spl', 0, ...
                    'frequency', [], ...
                    'samplingFrequency', GlobalSamplingFrequency, ...
                    'tspan', [Time(0), tf_extra]);
        case 'step'
            
            tf_extra = tf_extra + args.zeroDuration + args.fadeDuration;

            signal = ...
                StepSignal( ...
                    'amplitude', args.amplitude, ...
                    'frequency', args.frequency, ... defining characteristic position
                    'samplingFrequency', GlobalSamplingFrequency, ...
                    'zeroDuration', args.zeroDuration, ...
                    'fadeDuration', args.fadeDuration, ...
                    ...'default_analysis_start_time', args.zeroDuration, ...
                    'tspan', [Time(0), tf_extra]);
        case 'chirp'
            signal = Chirp( ...
                    'frequency', args.frequency, ...
                    args.amplitude_unit, args.amplitude, ...
                    'waveform', args.waveform, ...
                    'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
                    'onsetDuration', onsetDuration, ...
                    'offsetDuration', offsetDuration, ...
                    'zeroDuration', args.zeroDuration, ...
                    'fadeDuration', args.fadeDuration, ...
                    'tspan', [args.t0, signal_duration(tf_extra)]);
                    ...'tspan', [args.t0, tf_extra]);
        case 'greasy'
            [s, fs] = greasyWrapper(tf_extra);
            signal = VarSignal( ...
                'originalAudio', s, ...
                'originalSamplingFrequency', Frequency(fs, 'Hz'), ...
                'spl', args.amplitude, ...          [dB SPL]
                'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
                'tspan', [args.t0, tf_extra], ... Time(367.5, 'ms') % greasy
                'name', 'greasy');
        case 'click'
            signal = Click( ...
                args.amplitude_unit, args.amplitude, ...
                'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
                'onsetDuration', Time(1, 'us'), ...
                'offsetDuration', Time(1, 'us'), ...
                'clickDuration', Time(1e-3, 'us'), ...
                'zeroDuration', Time(1, 'ms'), ...
                'tspan', [args.t0, tf_extra]);
        case 'click2'
            signal = Click2( ...
                args.amplitude_unit, args.amplitude, ...
                'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
                'onsetDuration', Time(10, 'us'), ...
                'offsetDuration', Time(10, 'us'), ...
                'clickDuration', Time(1, 'us'), ...
                'zeroDuration', Time(1, 'ms'), ...
                'tspan', [args.t0, tf_extra]);
        case 'file'
            [s, fs] = audioread(args.signal_file);

            tf_extra = Time(numel(s) / fs, 's');

            [~, name, ~] = fileparts(args.signal_file);

            signal = VarSignal( ...
                'originalAudio', s, ...
                'originalSamplingFrequency', Frequency(fs, 'Hz'), ...
                'spl', args.amplitude, ...          [dB SPL]
                'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
                'onsetDuration', onsetDuration, ...
                'offsetDuration', offsetDuration, ...
                'tspan', [args.t0, tf_extra], ... Time(367.5, 'ms') % greasy
                'name', name);
        case 'var'
            s = args.signal_var.audio;
            fs = args.signal_var.fs;
            name = args.signal_var.name;

            assert(~isempty(s), 'Signal audio is empty')
            assert(~isempty(fs), 'Signal fs is empty')
            assert(~isempty(name), 'Signal name is empty')

            tf_extra = numel(s) / fs;

            signal = VarSignal( ...
                'originalAudio', s, ...
                'originalSamplingFrequency', fs, ...
                'spl', args.amplitude, ...          [dB SPL]
                'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
                'onsetDuration', onsetDuration, ...
                'offsetDuration', offsetDuration, ...
                'zeroDuration', args.zeroDuration, ...
                'fadeDuration', args.fadeDuration, ...
                'tspan', [args.t0, signal_duration(tf_extra)], ...
                'name', name);
        otherwise
            error('Unknown value %', args.special)
    end

    TOPT = topt_fun(tf_extra);

    if args.plotflag
        figure();
        plot(signal)
    end

    return

end

%%

signal = @(TSPAN) ...
        PureTone( ...
            'frequency', args.frequency, ...
            args.amplitude_unit, args.amplitude, ...
            'waveform', args.waveform, ...
            'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
            'onsetDuration', onsetDuration, ...
            'offsetDuration', offsetDuration, ...
            'zeroDuration', args.zeroDuration, ...
            'fadeDuration', args.fadeDuration, ...
            'tspan', TSPAN);
            % 'tspan', topt.total.tspan);

switch class(tf_extra)
    case 'Time' % tf_extra is the time                        
    case 'double' % tf_extra is the number of periods to simulate
        tf_extra = tf_extra / args.frequency;
    otherwise
        error('tf_extra must be either Time or double')
end

if ~args.strict_signal_length
    min_num_period = 4;
    if tf_extra * args.frequency < min_num_period
        % warning('Time of simulation seems too small');
    end
    tf_extra = max(tf_extra, min_num_period / args.frequency);
    
    min_num_time = Time(10, 'ms');
    tf_extra = max(tf_extra, min_num_time);
end

if ~args.strict_signal_length
    old_tf_extra = tf_extra;
    while true
        fs = GlobalSamplingFrequency.Hz;
        f = args.frequency.Hz;
        signal_length = tf_extra.s * fs;    
        N = signalLengthCheckSpectrumResolution(f, fs, signal_length);
        if N >= 5
            if old_tf_extra ~= tf_extra
    
                tf_extra = Time(ceil(tf_extra.ms), 'ms');
    
                if args.debug == true
                    dbprintf('Increasing tf_extra from %s to %s because of minimal spectrum resolution\n', ...
                        old_tf_extra, tf_extra)
                end
            end
            break
        else
            tf_extra = 1.1 * tf_extra;
        end
    end
end

% create time opt & stimulus
topt = topt_fun(tf_extra);
signal = signal(topt.total.tspan);

TOPT = topt;

% for i = 1:length(signal)
%     signal(i).eval();
%     figure
%     plot(signal(i))
% %     play(signal(i),'y')
% end
% return


end

function [ val ] = onset_period(onset, varargin)

if isa(onset, 'Time')
    val = onset;
    return
elseif isa(onset, 'char') || isa(onset, 'string')
    name = onset;
else
    error('onset must be of class Time, char or string, not %s', class(onset))
end

switch name
    case 'quick'
        val = quick_onset(varargin{:});
    case 'easy'
        val = easy_onset(varargin{:});
    case 'v1'
        val = onset_v1(varargin{:});
    otherwise
        error('not set up')
end
end

function [ onset ] = quick_onset(ampl, freq)
    arguments
        ampl (1,1) double
        freq (1,1) Frequency
    end
    onset_min_num_period = 2;
    onset = onset_min_num_period / freq;
end

function [ onset ] = easy_onset(ampl, freq)
    arguments
        ampl (1,1) double
        freq (1,1) Frequency
    end
    onset_min_num_period = 8 + 2.0*ampl;
    onset = onset_min_num_period / freq;
end

function [ onset ] = onset_v1(ampl, freq)
    arguments
        ampl (1,1) double
        freq (1,1) Frequency
    end
    onset = Time(max([10, ampl/5]), 'ms');
end