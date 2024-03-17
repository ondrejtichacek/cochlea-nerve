function [ TF ] = isequaleps( TOLERANCE, varargin )
%ISEQUALEPS

assert(numel(varargin) >= 2);

TF = true;

for i = 1:numel(varargin)-1
    A = varargin{i};
    B = varargin{i+1};
    
    TF = TF && all(abs(A - B) < TOLERANCE);
end

end

