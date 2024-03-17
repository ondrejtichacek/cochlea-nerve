function [ ] = optConsistency( PART, stimulus, mechopt, mnaopt, runopt, opt, plotopt, memopt, paropt )
%OPTCONSISTENCY 

switch PART
    case 'MECH'
    case 'MNA'
    case 'ANT'

        if mechopt.tf <= runopt.analysisStart
            warning('Total simulation time %f is too small for analysis to start (at %f). Consider extending the simulation', mechopt.tf, runopt.analysisStart);
        elseif mechopt.tf - runopt.analysisStart < runopt.analysisStart
            warning('Total simulation time %f is too small for long enough analysis starting at %f. Consider extending the simulation', mechopt.tf, runopt.analysisStart);
        end

    otherwise
end



end

