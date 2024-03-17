function [ N ] = signalLengthCheckSpectrumResolution(f, fs, signal_length, DELTA)
%SIGNALLENGTHCHECKSPECTRUMRESOLUTION 
arguments
    f
    fs
    signal_length
    DELTA = 0.02;
end

FAC = max([0.1,(1-log10(f/1000))]);

f_plus = exp(log(f)*(1+DELTA*FAC));
f_minus = exp(log(f)*(1-DELTA*FAC));

F = powerSpectrumResolution(signal_length, fs);
% xlim([f_minus, f_plus]);

f_plus_ind = find(F >= f_plus, 1);
f_minus_ind = find(F >= f_minus, 1);

assert(f_plus_ind > f_minus_ind)

N = (f_plus_ind - f_minus_ind);

end

