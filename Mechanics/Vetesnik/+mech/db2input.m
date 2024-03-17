function input = db2input(dB, AMo, L)
%DB2INPUT Rescales dB SPL into corresponding input values for human cochlea
%
% Original value was    AMo = 3e-6,     L = 0
% Now we use            AMo = 1e-6,     L = -10
%
% The parameter L could be incorporated in AMo or the other way round, but
% this way, it is more intuitive to calibrate the I/O curves.
%
% The parameter L essentailly shifts the I/O curve right (for L negative)
% or left (for L positive) by L db SPL.

input = AMo*10.^((dB + L)/20);

end