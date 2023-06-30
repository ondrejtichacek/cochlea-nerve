function [fun] = envelope_tolerance(stimulus, tspan, time_unit, fs, args)
%ENVELOPE_TOLERANCE
arguments
    stimulus 
    tspan (1,:) double
    time_unit (1,:) char
    fs (1,1) Frequency
    args.amplifier (1,1) logical = true
    args.plotflag (1,1) logical = false
    args.saturate_spl (1,1) double = +inf
end

if isa(stimulus, 'PureTone')

    if numel(stimulus) == 1
        envelope = stimulus.envelope_fcn(time_unit);
    else
        error('Not implemented')
    end

elseif isa(stimulus, 'Click') || isa(stimulus, 'Click2')
    envelope = stimulus.envelope_fcn_tolerance(time_unit);

elseif isa(stimulus, 'ZeroSignal')
    envelope = @(t) ones(size(t));

elseif isa(stimulus, 'CompoundSignal')
    envelope = stimulus.envelope_fcn(time_unit);
else
    envelope = @(t) ones(size(t));
    warning('Implement envelope')
end

%%

spl = stimulus.amplitude;
spl = max(spl); % TODO: this is not exact

fac = 1;
% fac = 10^(spl/20);
% if args.saturate_spl < +inf
%     fac = 10^(min([args.saturate_spl, spl])/20);
% end
fac = max(1, fac);

envelope_spl = @(t) spl + 20*log10(eps + abs(envelope(t)));
% envelope_spl = @(t) max(0, envelope_spl(t));

% gain_est_v0 = @(spl) max(0, 40 - 40*spl/60);

if args.amplifier == true
    g0 = 50;
else
    g0 = 0;
end

% g1 = 20;
gain_est_v1 = @(spl) boltzmann(-spl, g0, -50, 10);
% gain_est_v1 = @(spl) 0;
gain_est = gain_est_v1;

% if args.saturate_spl < +inf
%     gain_est = @(spl) gain_est(spl) - max([0, spl - args.saturate_spl]);
% end


% fac0 = 2;
% fac0 = 1;
% env = @(t) 10.^((-spl + envelope_spl(t) + gain_est(envelope_spl(t)))/20) / fac0;
env = @(t) envelope(t) .* (1 + 10.^(gain_est(envelope_spl(t))/20));
% env2 = @(t) max(10^(g1/20)/fac, env(t));

env2 = @(t) max(0.1, env(t));

figure
plot(env(tspan))

if args.plotflag
    h2 = figure;
    %tiledlayout flow
    %nexttile
    hold on
    s = linspace(0, 120, 200);
    %plot(s, gain_est_v0(s), 'DisplayName', 'v0')
    plot(s, gain_est_v1(s), 'DisplayName', 'v1')
    xlim([s(1), s(end)])
    xlabel('level (dB SPL)')
    ylabel('gain (dB)')
    legend()
    %nexttile
    %hold on
    % plot(s, s + gain_est_v0(s), 'DisplayName', 'v0')
    %plot(s, s + gain_est_v1(s), 'DisplayName', 'v1')
    %plot(s, s, 'Color', [0.65,0.65,0.65], 'DisplayName', 'passive')
    %xlim([s(1), s(end)])
    %ylim([s(1), s(end)])
    %xlabel('level (dB SPL)')
    %ylabel('response (dB)')
    %legend('Location','southeast')

    h1 = figure;
    tiledlayout(5,1)
    nexttile
    hold on
    plot(tspan*1e3, envelope(tspan), ...
        'DisplayName', 'envelope')
    plot(stimulus.audiotime.ms, stimulus.audio, ...
        'Color', [0.65,0.65,0.65], ...
        'DisplayName', 'waveform')
    ylabel('rel. amplitude (1)')
    legend()
    nexttile([2 1]);
    hold on
    plot(tspan*1e3, envelope(tspan), ...
        'DisplayName', 'envelope')
    plot(tspan*1e3, env(tspan), ...
        'DisplayName', 'envelope with gain')
    plot(tspan*1e3, env2(tspan), ...
        'DisplayName', 'tol_0')
    ylabel('rel. amplitude (1)')
    legend('Location','best')
    % xlabel('ms')
end

if isa(stimulus, 'CompoundSignal') || isa(stimulus, 'VarSignal')
    Ts = Time(12, 'ms'); % TODO: Improve
else
    Ts = (stimulus.Duration ...
        - stimulus.onsetDuration ...
        - stimulus.zeroDuration ...
        - stimulus.offsetDuration ...
        - stimulus.fadeDuration);
end
% convolution window
Tw = min(Ts, Time(24, 'ms'));

if isa(stimulus, 'PureTone')
    if Tw < Time(8, 'ms')
        error('Too short signal')
    end
end
if isa(stimulus, 'Click') || isa(stimulus, 'Click2')
    if Tw < Time(8, 'ms')
        Tw = Time(8, 'ms');
    end
end

ns = ceil(fs * Tw);

% gaussian
tw = linspace(0, Tw.s, ns);
w = gaussian(tw, 1, 0, Tw.s / 2);
w = [zeros(1,ns-101), ones(1,100), w];
% w = [zeros(1,ns-1), w];

% custom
% w = [zeros(1,ns-1), ns:-1:1];
% w = sqrt(w);
    
w = w ./ sum(w);

f = conv(env2(tspan), w, 'same');
ff = max(env2(tspan), f);

i = find(tspan > 8e-3 & env2(tspan) < 0.1, 1);
fff = zeros(size(tspan));
fff(i:end) = 0.1;
ff = max(fff, ff);

F = griddedInterpolant(tspan, ff);

if args.plotflag
    tt = linspace(tspan(1), tspan(end), 20);
    figure(h1);
    nexttile([2 1]);
    hold on    
    plot(tspan*1e3, f, '-', ...
        'DisplayName', 'tol_0 * window')
    plot(tspan*1e3, ff, '-', 'LineWidth', 2, ...
        'DisplayName', 'tolerance')
    plot(tspan*1e3, env2(tspan), '-', ...
        'DisplayName', 'tol_0 ')
    plot(tt*1e3, F(tt), 'ko', ...
        'DisplayName', 'interp')
    xlabel('time (ms)')
    ylabel('rel. amplitude (1)')
    legend('Location','best')

    axes('Position',[.78 .26 .1 .1])
    box on
    plot(1e3*[-flip(tw(2:end)),tw], w, 'k')
    ylim([-0.01, 0.115])
    yticks([])
    title('conv. w.')
    % xlabel('ms')

end

% v0
% fun = @(t) max(1/fac, envelope(t));
% v1
% fun = @(t) max(1/fac, env(t));
% v2
fun = @(t) F(t);

end

