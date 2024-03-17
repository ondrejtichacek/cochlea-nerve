function [ rgb ] = rwb( m )
%RWB
if nargin < 1, m = size(get(gcf,'colormap'),1); end

rgb = [ 202,0,32
        244,165,130
        94,60,153
        146,197,222
        5,113,176 ]/255;

% rgb = [ 8,29,88
%         37,52,148
%         34,94,168
%         29,145,192
%         65,182,196
%         127,205,187
%         199,233,180
%         237,248,177 % 237,248,177
%         255,255,204
%         255,237,160
%         254,217,118
%         254,178,76
%         253,141,60
%         252,78,42
%         227,26,28
%         189,0,38
%         128,0,38 ]/255;

if m < 4
    if m == 1
        rgb = rgb(3,:);
    elseif m == 2
        rgb = rgb([1,5],:);
    elseif m == 3
        rgb = rgb([1,3,5],:);
    end
else
    rgb = interp1(linspace(0,1, size(rgb,1)),rgb(:,:),linspace(0,1,m));
end

% b = repmat(linspace(0,1,200),20,1);
% imshow(b,[],'InitialMagnification','fit')

end

