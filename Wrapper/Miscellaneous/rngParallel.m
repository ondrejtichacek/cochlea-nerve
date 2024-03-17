function [ ] = rngParallel( RNG_worker, RNG_local, pool, seed )
%RNGPARALLEL 

% shuffle local RNG, seed RNG on workers, set worker RNG algorithm to Mersene-Twister
% Example rngParallel('twister','shuffle');

if nargin >= 2    
    rng(RNG_local);
end

if nargin < 3
    pool = gcp();
end

if nargin < 4
    seed = randi(floor(intmax/10),pool.NumWorkers,1);
end

% Shuffle on each parallel worker
spmd
    rng(seed(labindex),RNG_worker);
end


end

