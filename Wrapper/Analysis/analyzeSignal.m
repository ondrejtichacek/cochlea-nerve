function analyzeSignal( ...
    VAR, REF, x, filename, stimulus, samplingFrequency, mechopt, args)

arguments
    VAR  (:,:)  double
    REF  (1,:)  double
    x    (1,:)  double
    filename    char
    stimulus    
    samplingFrequency Frequency
    mechopt     mechOpt
    args.do_powerSpectrum = true
    args.do_HilbertPhase = true
    args.do_std_width = true
    args.do_travelling_wave = true
    args.do_acdc = true
    args.do_oscillation_argmax_std = false
    args.do_oscillation_char_frequency_model = false
    args.do_oscillation_max_delta_restricted = false
    args.do_oscillation_spectrum_restricted = true
end

stats = struct();

tdim = 1;
ndim = 2;

n = size(VAR, ndim);
nt = size(VAR, tdim);


%% calculate power spectrum at each cross-section
if args.do_powerSpectrum == true

    VAR_REL = VAR - mean(VAR, 1);
    
    % find peaks in the spectrum to note main frequencies
    f = powerSpectrum(VAR_REL(:,1), samplingFrequency.Hz);
    P = zeros(numel(f), n);
    % F{n} = [];
    [max_freq, max_freq_power, max_freq_ind] = deal(zeros(1,n));
    for y = 1:n
        [~, P(:,y)] = powerSpectrum(VAR_REL(:,y), samplingFrequency.Hz);
        % [~,F{y}] = findpeaks(P(:,y),f,'SortStr','descend');
        [max_freq_power(y), max_freq_ind(y)] = max(P(:,y));
        max_freq(y) = f(max_freq_ind(y));        
    end

    stats.power_spectrum = struct( ...
        'P', P, ...
        'x', x, ...
        'max_freq', max_freq, ...
        'max_freq_ind', max_freq_ind, ...
        'max_freq_power', max_freq_power, ...
        'freq', f);
    
    % how to plot:
    % [X,Y] = meshgrid(f,x);
    % mesh(Y,X,P');
    % ylim([0,5000])
    
end

% Calculate the phase of the variable using Hilbert transfrom
if args.do_HilbertPhase == true
    % at every step, the spatial quantity is transformed
    % and finally averaged over all the time steps.
    [Signal_Toolbox_licenceAvailable, ~] = license('checkout','Signal_Toolbox');
    if Signal_Toolbox_licenceAvailable
        phase = zeros(1,n);
        t1 = 1;
        t2 = nt;
        for tstep = 1:(t2 - t1 + 1)
            Hs = hilbert(VAR(tstep,:));
            ph = unwrap(-angle(Hs));
            ph = ph - ph(1);
            phase = phase + ph;
        end
        phase = phase/(t2 - t1 + 1);
    else
        warning('Can''t use hilbert transform without signal processing toolbox. Not including phases in the ouput.')
        phase = {};
    end

    stats.phase = phase;
end

%% Save some useful information
stats.x = x;

ver = mechopt.char_frequency_model_best_fit;
[cx, cf] = characteristic_frequency_model(ver);

has_frequency = ~isempty(stimulus.frequency);

main_freq = stimulus.frequency;

if numel(stimulus.frequency) > 1
    if all(stimulus.frequency == stimulus.frequency(1))
        % all components have same frequency
        main_freq = stimulus.frequency(1);
    else
        warning('Analysis not implemented for multiple signals')
        has_frequency = false;
    end
end

if has_frequency
    stats.char_position = cx(main_freq.Hz);
    stats.char_position_ind = ceil(stats.char_position * n);
end

%% Simple statistics
stats.max = max(VAR(:));

S = std(VAR, 0, tdim);
stats.std = S;

if args.do_std_width == true
    [Signal_Toolbox_licenceAvailable, ~] = license('checkout','Signal_Toolbox');
    if Signal_Toolbox_licenceAvailable
        [~,~,stats.std_width,~] = findpeaks(S, x, 'Npeaks', 1, 'WidthReference', 'halfheight');
    end
end

stats.time_mean = mean(VAR, tdim);

%% Profile
if args.do_travelling_wave == true
    stats.TWupper = max(VAR, [], tdim);
    stats.TWlower = min(VAR, [], tdim);
end
if args.do_acdc == true
    stats.DC = abs(mean(VAR, tdim));
    stats.AC = max(abs(VAR - stats.DC), [], tdim);
end
%% Oscillation maxima

if has_frequency
%% Find oscillation maximum using power spectrum restricted

if args.do_oscillation_spectrum_restricted

    [~, ~, PP] = analysis.oscillation_spectrum_core( ...
            stats.power_spectrum.P, ...
            stats.power_spectrum.freq, ...
            main_freq.Hz);

    assert(size(PP,2) == mechopt.Numstacks);
    [~, ind] = max(max(PP, [], 1), [], 2);

    stats.oscillation_max.spectrum_restricted = oscillationStatistics(VAR, ind);
end

%% Find oscillation maximum using std
if args.do_oscillation_argmax_std == true
    [~, ind] = max(S);

    stats.oscillation_max.argmax_std = oscillationStatistics(VAR, ind);
end

%% Find oscillation maximum using characteristic frequency model
if args.do_oscillation_char_frequency_model == true
    if isa(stimulus, 'PureTone')

        % characteristic position
        % cf_model = 'loglinear';
        % cf_parameters = [4.3647; -2.8577];
        % cx = characteristic_position(stimulus.frequency.Hz, cf_model, cf_parameters);

        if mechopt.gain == 0 || (isnumeric(mechopt.lambda) && all(mechopt.lambda == 0))
            model = 'passive';
        else
            model = 'active';
        end

        error('first fix warning when loading char freq model');

        cx_model = characteristic_frequency_model(model, mechopt.middle_ear_params.identifier);
        cx = cx_model(stimulus.amplitude, main_freq.Hz);

        [~, ind] = find(x >= cx, 1);

        stats.oscillation_max.char_frequency_model = oscillationStatistics(VAR, ind);
    end
end

%% Find oscillation maximum using maximal delta restricted by characteristic frequency model
if args.do_oscillation_max_delta_restricted == true
    if isa(stimulus, 'PureTone')

        % The actual maximal position is in the interval
        %   [xx - 0.1, xx + 0.05]
        % where xx is the characterisitic position of the frequency map.
        % This is a heuristic result
        i1 = find(x >= stats.char_position - 0.1, 1);
        i2 = find(x >= stats.char_position + 0.05, 1);

        % heuristic metric to find maxima of oscillations: max(abs(V - V0)) + 2*std(V)
        [~, ind] = max(2 * sqrt(S(:, i1:i2)) + max(abs(VAR(:, i1:i2) - REF(:, i1:i2)), [], tdim) );

        ind = i1 + ind - 1;

        stats.oscillation_max.max_delta_restricted = oscillationStatistics(VAR, ind);
    end
end

%%
end

function s = oscillationStatistics(VAR, ind)
    s = struct();

    s.ind = ind;

    s.char_frequency = Frequency(cf(x(ind)), 'Hz');

    s.xpos = x(ind);
    s.val = VAR(:, ind);
    s.max = max(VAR(:,ind));
    s.min = min(VAR(:,ind));

    if args.do_powerSpectrum == true
        [s.main_freq, s.trust_this_data] = checkFrequency(ind);
    end
end

function [main_freq, trust_this_data] = checkFrequency(ind)
    % Compare the frequency of signal change at the maximal-variance
    % place with the frequency of the stimulus. If the stimulus is a
    % pure tone, they should approximately match.

    switch class(stimulus)
    case {'PureTone', 'CompoundSignal'}

        if isa(stimulus, 'PureTone')
            F = stimulus.frequency;
        elseif isa(stimulus, 'CompoundSignal') && all(stimulus.frequency == stimulus.frequency(1))
            F = stimulus.frequency(1);
        else
            error("Can't do that yet")
        end
        

        main_freq = Frequency(max_freq(ind), 'Hz'); % frequency with the highest power
        % main_freq = stats.char_frequency; % characteristic frequency from the frequency map
        
        rel_freq_difference = abs(main_freq.Hz - F.Hz)/F.Hz;
        
        if rel_freq_difference > 0.05
            warning(['Frequency of stimulation at the maximum-amplitude position was %s, ', ...
                    'stimulus frequency %s. Possible error (artificial oscillations)? ', ...
                    'ind = %d, rel pos = (%g)'], ...
                main_freq, F, ind, ind/n)
            trust_this_data = false;
        else
            trust_this_data = true;
        end
        
    case 'ZeroSignal'
        main_freq = nan;
        trust_this_data = true;
        
    otherwise
        error('We did not perform check for frequency of displacement because it might not make sense for the class %s. Only valid for pure-tone signals.', class(stimulus))
    end
end

%% save space

stats.power_spectrum.P = [];

%% Save results

if exist(filename, 'file')
    sdelete(filename);
end

save(filename, '-struct', 'stats', '-v7.3', '-nocompression')

end