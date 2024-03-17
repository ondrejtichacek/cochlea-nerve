function [ Y ] = movingAverageFilter( windowSize, SIGNAL )
%MOVINGAVERAGEFILTER
arguments
    windowSize (1,1) double
    SIGNAL (:,:) double
end

% moving average filter
% windowSize = 500;
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;        
% 
% Y = filter(b,a, [ SIGNAL(1) * ones(windowSize,1); SIGNAL(:) ]);
% Y = Y((windowSize+1):end);


% convolution based version

[nt, ns] = size(SIGNAL);

w = (1/windowSize)*ones(1,windowSize);

Y = zeros(nt, ns);

for i = 1:ns
    tmp = conv([ ...
            SIGNAL(1,i)*ones(1,2*windowSize), ...
            SIGNAL(:,i)', ...
            SIGNAL(end,i)*ones(1, 2*windowSize)], ...
        w, 'same');
    Y(:,i) = tmp((2*windowSize+1):end-(2*windowSize));
end

end

