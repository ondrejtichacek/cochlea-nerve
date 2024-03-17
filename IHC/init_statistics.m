function v = init_statistics(xind, n, sr, PeakFields, time, Signal_Toolbox_licenceAvailable)
    
    v.xind = xind;
    v.sr = sr;

    v.phase = [];
    
    v.N = 0;
    [ v.MEAN, ...
      v.VAR, ...
      v.SD, ...
      v.M2 ] = deal(zeros(numel(time), n));
  
    v.MIN = Inf(numel(time), n);
    v.MAX = -Inf(numel(time), n);

    if Signal_Toolbox_licenceAvailable == 1

        v.InterTimes = [];
        v.Peak = [];
        v.Location = [];
        v.Width = [];
        v.Prominence = [];

        for f = 1:length(PeakFields)

            v.(PeakFields{f}).N = 0;
            [ v.(PeakFields{f}).MEAN, ...
            v.(PeakFields{f}).M2 ] = deal(zeros(1,n));

            v.(PeakFields{f}).MIN = Inf(1,n);
            v.(PeakFields{f}).MAX = -Inf(1,n);

        end
    end
    
    OtherStatistics = {'TimeAverage'};
    for f = 1:length(OtherStatistics)

        v.(OtherStatistics{f}).N = 0;
        [ v.(OtherStatistics{f}).MEAN, ...
        v.(OtherStatistics{f}).M2 ] = deal(zeros(1,n));

        v.(OtherStatistics{f}).MIN = Inf(1,n);
        v.(OtherStatistics{f}).MAX = -Inf(1,n);

    end

    v.PHASE = [];
end
