function [sig] = band_noise_bursts(duration, fcentre, fs, gap_duration)
%BAND_NOISE_BURSTS 

% fc_lims = Frequency([6,16], 'kHz');
% bw_lims = Frequency([2,4], 'kHz');

num_bursts = numel(fcentre);
noise_seeds = int32(linspace(0, double(intmax), num_bursts));

sig = [];

rng_state = rng();

for i = 1:num_bursts

    % fc = linf(fc_lims(1), fc_lims(2), rand(1));
    % bw = linf(bw_lims(1), bw_lims(2), rand(1));

    fd = 2^(1/6); % third octave bands
    fupper = fcentre(i) * fd;
    flower = fcentre(i) / fd;

    fc = mean([flower, fupper]);
    bw = fupper - flower;

    rng(noise_seeds(i))

    s = sig_bandpassnoise(fc.Hz, fs.Hz, duration.s, 0, bw.Hz);
    
    s = s(:)';

    ns = numel(s);
    t = Time(linspace(0,ns * 1/fs.Hz, ns), 's');
    
    t0 = Time(0, 'ms');
    t3 = duration;
    e1 = Time(2, 'ms');
    e2 = Time(2, 'ms');
    env = Signal.envelope_cos_core(t.s, t0.s, t3.s, e1.s, e2.s);

    sig = [sig, s .* env];

    gap = zeros(1, gap_duration * fs);

    if i < num_bursts
        sig = [sig, gap];
    end

end

rng(rng_state)

end

