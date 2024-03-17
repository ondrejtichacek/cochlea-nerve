function [ AnalyticSignal, SignalAmplitude, SignalPhase, InstantaneousFrequency ] = HilbertTransform( Signal, Fs )
%HILBERTTRANSFORM

% hibert returns "analytic signal" (complex)
AnalyticSignal = hilbert(Signal);

% local amplitude information (envelope)
SignalAmplitude = abs(AnalyticSignal);

% phase information (wrapped on unit circle)
SignalPhase = angle(AnalyticSignal);

% instantaneous frequency
InstantaneousFrequency = diff(unwrap(SignalPhase))/((1/Fs.Hz)*2*pi);

end

