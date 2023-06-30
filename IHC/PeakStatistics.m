function [ STAT, InterTimes ] = PeakStatistics( TIME, pks, locs, width, prom )
%PEAKSTATISTICS 
arguments
    TIME
    pks
    locs
    width = []
    prom = []
end

if isempty(TIME) && isempty(pks) && isempty(locs)
    
    STAT = { ...
            'NumberOfPeaks', ...
            'MinPeak', 'MaxPeak', 'MeanPeak', 'StdPeak', ...
            'MinProminence', 'MaxProminence', 'MeanProminence', 'StdProminence', ...
            'MinWidth', 'MaxWidth', 'MeanWidth', 'StdWidth', ...
            'MinInterval', 'MaxInterval', 'MeanInterval', 'StdInterval', ...
            'MinPeak', 'MaxPeak', 'MeanPeak', 'StdPeak', ...
            'GlobalPeakRate' ...
    };

    return

end

STAT.NumberOfPeaks = numel(pks);

% peak value
STAT.MinPeak     = min(pks);
STAT.MaxPeak     = max(pks);
STAT.MeanPeak    = mean(pks);
STAT.StdPeak     = std(pks);

% peak prominence
if ~isempty(prom)
    STAT.MinProminence   = min(prom);
    STAT.MaxProminence   = max(prom);
    STAT.MeanProminence  = mean(prom);
    STAT.StdProminence   = std(prom);
end

% peak width (half-prominence)
if ~isempty(width)
    STAT.MinWidth    = min(width);
    STAT.MaxWidth    = max(width);
    STAT.MeanWidth   = mean(width);
    STAT.StdWidth    = std(width);
end

% inter-peak interval
if length(locs) > 1
    InterTimes = locs(2:end) - locs(1:end-1);
else
    InterTimes = zeros(0,1);
end

STAT.MinInterval    = min(InterTimes);
STAT.MaxInterval    = max(InterTimes);
STAT.MeanInterval   = mean(InterTimes);
STAT.StdInterval    = std(InterTimes);

% peak global rate
STAT.GlobalPeakRate = STAT.NumberOfPeaks / (TIME(end) - TIME(1)) *1e3;

end

