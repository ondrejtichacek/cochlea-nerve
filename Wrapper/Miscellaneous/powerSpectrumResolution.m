function [ F ] = powerSpectrumResolution( signal_length, fs )
arguments
    signal_length
    fs
end

m = signal_length;       % original sample length
n = pow2(nextpow2(m));   % transform length
f = (0:n-1)*(fs/n);
F = f(1:floor(n/2));

end

