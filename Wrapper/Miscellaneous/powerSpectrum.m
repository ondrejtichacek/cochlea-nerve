function [ F, P ] = powerSpectrum( signal, fs, args )
arguments
    signal
    fs
    args.plotflag (1,1) logical = false
end

m = length(signal);       % original sample length
n = pow2(nextpow2(m));    % transform length
y = fft(signal,n);        % DFT of signal

f = (0:n-1)*(fs/n);
power = abs(y).^2/n;

F = f(1:floor(n/2));
P = power(1:floor(n/2));

if args.plotflag == true
    figure;
    plot(F, P);
end

end

