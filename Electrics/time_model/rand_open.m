function [open] = rand_open(probability, N)
%RAND_OPEN
arguments
    probability (:,1) double
    N (1,1) double {mustBeInteger}
end

K = numel(probability);

open = sum(rand(K, N) < probability, 2) / N;

end

