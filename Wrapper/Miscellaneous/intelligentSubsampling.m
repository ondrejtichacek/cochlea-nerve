function [ S, T, IND ] = intelligentSubsampling( SIGNAL, TIME, FACTOR, varargin )
%INTELIGENTSUBSAMPLING

%% V1 (much simpler, does not support all varargin)

% [~,locsMAX] = findpeaks(SIGNAL, varargin{:});
% [~,locsMIN] = findpeaks(-SIGNAL, varargin{:});
% 
% IND = sort(unique([1:FACTOR:length(TIME), length(TIME), locsMAX, locsMIN]));
% 
% S = SIGNAL(IND);
% T = TIME(IND);

%% V2

for i = 1:numel(varargin)
    if isa(varargin{i}, 'Unit')
        varargin{i} = varargin{i}.(varargin{i}.si_unit);
    end
end

if isa(TIME, 'Unit')
    t = TIME.(TIME.si_unit);
else
    t = TIME;
end

if isa(SIGNAL, 'Unit')
    s = SIGNAL.(SIGNAL.si_unit);
else
    s = SIGNAL;
end

[~,locsMAX] = findpeaks(s, t, varargin{:});
[~,locsMIN] = findpeaks(-s, t, varargin{:});

L = sort(unique([locsMAX(:); locsMIN(:)]));

locsIND = zeros(1,length(L));
for i = 1:length(L)
    locsIND(i) = find(L(i) == t);
end

IND = sort(unique([1:FACTOR:length(t), length(t), locsIND]));

S = SIGNAL(IND);
T = TIME(IND);

end

