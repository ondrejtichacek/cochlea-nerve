function [ stat ] = ReplicationStatistics( x, Time, stat, PEAK_PROPERTIES, Signal_Toolbox_licenceAvailable, args )
arguments
    x
    Time
    stat
    PEAK_PROPERTIES
    Signal_Toolbox_licenceAvailable

    args.do_find_peaks (1,1) logical = true
end

[stat.N, stat.MEAN, stat.M2, stat.VAR, stat.SD] = OnlineVariance( x, stat.N, stat.MEAN, stat.M2 );

[stat.MIN, stat.MAX] = OnlineMinMax( x, stat.MIN, stat.MAX );

if args.do_find_peaks && Signal_Toolbox_licenceAvailable == 1
    
    warning('error','signal:findpeaks:largeMinPeakHeight'); % http://undocumentedmatlab.com/blog/trapping-warnings-efficiently
    try
        [pks,locs,width,prom] = findpeaks( x, Time, PEAK_PROPERTIES{:});
    catch ME
        switch ME.identifier
            case 'signal:findpeaks:largeMinPeakHeight'
                [pks,locs,width,prom] = deal(zeros(0,1));
            otherwise
                throw(ME);
        end
    end
    warning('on','signal:findpeaks:largeMinPeakHeight');
    
    [STAT, InterTimes] = PeakStatistics( Time, pks, locs, width, prom );

    PeakFields = fields(STAT);

    stat.InterTimes{end+1} = InterTimes;
    stat.Peak{end+1} = pks;
    stat.Location{end+1} = locs;
    stat.Width{end+1} = width;
    stat.Prominence{end+1} = prom;


    for f = 1:length(PeakFields)

        [ stat.(PeakFields{f}).N, ...
          stat.(PeakFields{f}).MEAN, ...
          stat.(PeakFields{f}).M2, ...
          stat.(PeakFields{f}).VAR, ...
          stat.(PeakFields{f}).SD ] = OnlineVariance( STAT.(PeakFields{f}), ...
                                                   stat.(PeakFields{f}).N, ...
                                                   stat.(PeakFields{f}).MEAN, ...
                                                   stat.(PeakFields{f}).M2 );
        [ stat.(PeakFields{f}).MIN, ...
          stat.(PeakFields{f}).MAX ] = OnlineMinMax( STAT.(PeakFields{f}), ...
                                                   stat.(PeakFields{f}).MIN, ...
                                                   stat.(PeakFields{f}).MAX );
    end
    
%     figure;
%     hold on
%     plot(Time, x)
%     plot(locs, pks, 'x')
    
    
end

OtherStatistics = struct( ...
    'TimeAverage', mean(x, 1));

PeakFields = fields(OtherStatistics);


for f = 1:length(PeakFields)

    [ stat.(PeakFields{f}).N, ...
      stat.(PeakFields{f}).MEAN, ...
      stat.(PeakFields{f}).M2, ...
      stat.(PeakFields{f}).VAR, ...
      stat.(PeakFields{f}).SD ] = OnlineVariance( OtherStatistics.(PeakFields{f}), ...
                                               stat.(PeakFields{f}).N, ...
                                               stat.(PeakFields{f}).MEAN, ...
                                               stat.(PeakFields{f}).M2 );
    [ stat.(PeakFields{f}).MIN, ...
      stat.(PeakFields{f}).MAX ] = OnlineMinMax( OtherStatistics.(PeakFields{f}), ...
                                               stat.(PeakFields{f}).MIN, ...
                                               stat.(PeakFields{f}).MAX );
end
   
    

end
