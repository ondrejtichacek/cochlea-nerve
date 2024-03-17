function [ N, MEAN_, M2, VAR_, SD_ ] = OnlineVariance( x, N, MEAN_, M2 )
%ONLINEVARIANCE
%
% Calculates mean, var and std for online data, i.e. does not require all
% data in memmory at once.
%
% USAGE:
% assume data:
%
% data = rand(2,4,8,16);
%
% then to calculate
%
% mean(data,3)
% var(data,0,3)
% std(data,0,3)
%
% it is enough to do the following
%
% [N, MEAN, M2] = deal([]);
% for i = 1:size(data,3)
%     [ N, MEAN, M2, VAR, SD ] = OnlineVariance( data(:,:,i,:), N, MEAN, M2 );    
% end
%
%
%
% TEST SCRIPT
% see OnlineStatsTest

if nargin == 1 || isempty(N) || isempty(MEAN_) || isempty(M2)
    N = 0;
    MEAN_ = zeros(size(x));
    M2 = zeros(size(x));    
end

if isempty(x) || all(isnan(x(:)))
    % return inputs        
else    
    if any(isnan(x))
        error('Don''t know how to process arrays containing some (and not all) NaNs');
    else
        N = N + 1;

        delta = x - MEAN_;
        MEAN_ = MEAN_ + delta/N;
        M2 = M2 + delta .* (x - MEAN_);
            
    end
end

if nargout > 3
    if N >= 2
        VAR_ = M2 / (N - 1);
        SD_ = sqrt(VAR_);
    else
        VAR_ = NaN;
        SD_ = NaN;
    end
end

end

