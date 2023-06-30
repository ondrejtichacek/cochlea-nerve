function [status] = MEplotting( resFiles, plotopt, stimulus, topt, midopt, runopt, opt, memopt )
%MEPLOTTING
arguments
    resFiles (1,1) struct
    plotopt (1,1) plotOpt
    stimulus
    topt (1,1)
    midopt (1,1) midOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
end

doplot = plotopt.do.middle_ear;

%%

status = true;

% matfile alias

if doplot.output == true

    fY = resFiles.middle_ear.fy;
    fT = resFiles.middle_ear.ft;

    t = Time(fT.t, midopt.simulation_units.time);
    Y = fY.y(:,midopt.IND);

    hfig = figure();
    plot(t.ms, Y)
    xlabel('time (ms)')

end



end

