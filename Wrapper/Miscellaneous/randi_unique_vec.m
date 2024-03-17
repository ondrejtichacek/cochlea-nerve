function [X] = randi_unique_vec(imax, n)

imin = 1;

if numel(imax) == 2
    imin = imax(1);
    imax = imax(2);
end

assert(imax > imin, 'imax must be > imin')
assert(imax - imin + 1 >= n, 'imax - imin must be >= n')

X = imin - 1 + randperm(imax - imin + 1, n);

% X = [];
% while numel(unique(X)) ~= n
%     X = randi([imin, imax], n, 1);
% end

end

