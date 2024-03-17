function [f, SS, SS_norm] = plotSignalFT( y, Fs )
%PLOTSIGNALFT 
% y  ... signal
% Fs ... sampling frequency

L = length(y);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y

S = fft(y,NFFT)/L;
SS = 2*abs(S(1:NFFT/2+1));

y_norm = y / max(abs(y));

S_norm = fft(y_norm,NFFT)/L;
SS_norm = 2*abs(S_norm(1:NFFT/2+1));


f = Fs/2*linspace(0,1,NFFT/2+1);

end

