function [cmp] = lcmp(N)
%LCMP 

if N == 1
    cmp = [0,0,0];
else
    cmp = parula(2*N);
    cmp = cmp(N-1:end-2,:);
    cmp = flipud(cmp);
end
end

