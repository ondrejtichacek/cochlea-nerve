clear all

externalResonanceFilters = [ ... gain order lp hp
  -11.99927987142523,  2,   206.2119974669669,   358.5102170379319;
  -11.99998556203806,  2,   288.0827584566169,   444.8320921474029;
  -6.660289825314575,  2,    604.601203475762,   1656.284925736281;
  -6.009994624924989,  2,   365.6631053098679,   1177.111851175218;
  -4.392649142756406,  2,   1600.417554322875,   3593.838940335973;
  -3.819930936639489,  2,   1634.013499558196,   3583.542619484177;
 -0.9495775823395274,  2,   2415.056217002997,   5319.282950031623;
   8.077026047557176,  2,   2435.063884295232,   4541.144792545283;
   3.016053335238496,  2,   2809.307899606387,    7618.45593491728;
   6.209711798089334,  2,   1398.641611236509,   4696.586038747771;
  -9.118556974337517,  2,   6461.689628623291,   9721.286323027263;
];


% FREQUENCY = round(erbspace(100,20000,64));
FREQUENCY = round(erbspace(100,20000,256));

max_pressure = zeros(size(FREQUENCY));
phase_lag_hilbert = zeros(size(FREQUENCY));

% fs = Frequency(2000, 'kHz');
fs = Frequency(400, 'kHz');

sampleRate = fs.Hz;

dt = 1/fs;
dt = dt.s;

t = 0:dt:0.5;

figure
hold on
grid on
set(gca, 'XScale', 'log')
erbticks;
xlim([125,20000]);
ylim([-100, 15]);
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')

for i = 1:numel(FREQUENCY)

    f = Frequency(FREQUENCY(i), 'Hz');

    level = 0;

    p0 = 1;
    p = p0 * 10.^(level/20);

    inputPressure = p * sin(2*pi*f.Hz*t);

    ExtEarPressure = OuterEar_Meddis(inputPressure, sampleRate, ...
        'externalResonanceFilters', externalResonanceFilters);

    n = numel(ExtEarPressure);
    nn = ceil(n/2);
    
    inputPressure = inputPressure(nn:end);
    ExtEarPressure = ExtEarPressure(nn:end);
    
    max_pressure(i) = max(ExtEarPressure);
    
    phase_lag_hilbert(i) = mean( ...
          unwrap(angle(hilbert(ExtEarPressure))) ...
        - unwrap(angle(hilbert(inputPressure))));
    
    if mod(i,16) == 0
        fres = max(20, f.Hz / 20);
        [pxx,f] = pspectrum(ExtEarPressure, sampleRate, ...
            'FrequencyResolution', fres);

        plot(f,pow2db(pxx))

        [pxx,f] = pspectrum(inputPressure, sampleRate, ...
            'FrequencyResolution', fres);

        plot(f,pow2db(pxx))

        p0;
    end
    
end

magnitude = 20*log10(max_pressure/max(inputPressure));

%%

runopt.save_figures = false;
% runopt.save_figures = true;

hfig = figure;

hold on
grid on

set(gca, 'XScale', 'log')
erbticks;

xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

[x, ya, yb] = Mehrgardt_1977.fig20();

plot(x, yb);
plot(x, ya);

plot(FREQUENCY, magnitude, 'k');

legend({'Shaw (1974)', 'Mehrgardt (1977)', 'model'}, ...
    'Location', 'southwest')

% Save figure

FIGNAME = sprintf('OE-magnitude-dev');
PATH = 'img/OE/';

export_flag = false;
% export_flag = true;

if runopt.save_figures
    analysis.plot.export_figure(hfig, FIGNAME, PATH, export_flag, 'Results/OE')
end

%%


hfig = figure;

hold on
grid on

set(gca, 'XScale', 'log')
erbticks;

xlabel('Frequency (Hz)')
ylabel('Phase (deg)')

plot(FREQUENCY, 180/pi*unwrap(phase_lag_hilbert), 'k');

% legend({'model'})

% Save figure

FIGNAME = sprintf('OE-phase-dev');
PATH = 'img/OE/';

export_flag = false;
% export_flag = true;

if runopt.save_figures
    analysis.plot.export_figure(hfig, FIGNAME, PATH, export_flag, 'Results/OE')
end