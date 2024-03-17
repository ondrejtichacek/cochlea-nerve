function [ V, T, S ] = ANTPreprocess( Voltage, TimeVoltage, Signal, TimeSignal, antopt, topt )
%ANTPREPROCESS

Signal = Signal(:);
TimeSignal = TimeSignal(:);

[ntV, nX] = size(Voltage);
[ntS, ~] = size(Signal);

assert(all(size(TimeSignal) == [ntS,1]), 'TimeSignal size is not valid');
assert(all(size(TimeVoltage) == [ntV,1]), 'TimeVoltage size is not valid');

%% TMP fix

Signal(isnan(Signal)) = 0;

%% INTERPOLATE

t_eps = Time(1, 'ns');

if true
    % t0 = TimeVoltage(1);
    % tf = TimeVoltage(end);
    
    t0 = topt.total.t0;
    tf = topt.total.tf;
    
    assert(t0 >= TimeVoltage(1) - t_eps, ...
        ['Requested time points are not consistent with Voltage Time from MNA;\n', ...
         'debug info: (t0 - TimeVoltage(1) = %s)'], ...
        t0 - TimeVoltage(1))
    assert(tf <= TimeVoltage(end) + t_eps, ...
        ['Requested time points are not consistent with Voltage Time from MNA;\n', ...
         'debug info: (tf - TimeVoltage(end) = %s)'], ...
        tf - TimeVoltage(end))
    
    signalLength = tf - t0; % [s]
    numberOfPoints = 1 + ceil(antopt.samplingFrequency * signalLength);
    T = linspace( t0, tf, numberOfPoints ); % time points
    T = T';
else
    T = TimeVoltage;
end

if nX == 1
    % V = interp1(TimeVoltage.us, Voltage, T.us, 'linear', 'extrap');
    F = griddedInterpolant(TimeVoltage.us, Voltage);
    V = F(T.us);
else
    F = griddedInterpolant({TimeVoltage.us, 1:nX}, Voltage);
    V = F({T.us, 1:nX});
end

% S = interp1(TimeSignal.us, Signal, T.us, 'linear', 'extrap');
F = griddedInterpolant(TimeSignal.us, Signal);
S = F(T.us);

assert(size(V,1) == length(T), 'Interpolated voltage has wrong shape')
assert(size(S,1) == length(T), 'Interpolated signal has wrong shape')

%% Force quick onset (for testing purposes only)

if antopt.force_quick_onset == true
    
    error('this does not work well for high frequency signals')

    warning('Onset of the IHC voltage essentially deleted.')
    
    assert(nX == 1);

    Vss = V(1,:);
    Vrel = V - Vss;
    A = max(Vrel,[], 1);
    i = find(Vrel > A*0.9, 1);
    j = find(Vrel(i:end,:) < 0, 1);
    ind = i+j;
    V(1:ind,:) = Vss;

end

%% SCALE

FACTOR = antopt.VoltageAmplitudeFactor;

if FACTOR ~= 1

%     m = mean(V);
%     m = -55.14;
    m = V(:,0);

    V = FACTOR*(V - m) + m;
    warning('Relative voltage amplitude modified by factor %g', FACTOR)
    
end

%% SCALE 2

FACTOR = antopt.VoltageAmplitudeAddFactor;

if FACTOR ~= 0

%     m = mean(V);
%     m = -55.14;
    m = V(:,0);

    M = max(V-m, [], 1);
    
    for i = 1:numel(M)
        
        V(:,i) = (M(i) + FACTOR) / M(i) * (V(:,i) - m) + m;
        
    end
    
    warning('Relative voltage amplitude modified by additive factor %g', FACTOR)
    
end

end

