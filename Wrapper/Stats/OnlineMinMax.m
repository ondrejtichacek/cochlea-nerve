function [ MIN_, MAX_ ] = OnlineMinMax(x, MIN_, MAX_)
%ONLINEMINMAX%
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
% min(data,[],3)
% max(data,[],3)
%
% it is enough to do the following
%
% [MIN, MAX] = deal([]);
% for i = 1:size(data,3)
%     [ MIN, MAX ] = OnlineMinMax( data(:,:,i,:), MIN, MAX );
% end
%
%
%
% TEST SCRIPT
% see OnlineStatsTest

if nargin == 1 || isempty(MIN_) || isempty(MAX_)
	MIN_ = Inf(size(x));
    MAX_ = -Inf(size(x));
end

sel_min = find(x < MIN_);
sel_max = find(x > MAX_);

MIN_( sel_min ) = x( sel_min );
MAX_( sel_max ) = x( sel_max );

end