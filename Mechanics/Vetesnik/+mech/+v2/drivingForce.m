function [ g ] = drivingForce( t, N, in, TMa, ~, stimulus_fcn, ~ )
%DRIVINGFORCE
% t, N, in, TMa, BMy, stimulus, YCut

inp = -stimulus_fcn(t);
g = [zeros(N,1); in*inp; zeros(N,1); -TMa.*in*inp];

end

