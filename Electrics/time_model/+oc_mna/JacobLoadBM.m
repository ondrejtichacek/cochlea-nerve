function [vars] = JacobLoadBM( tspan, fMechResult, unit, xgrid_mech, xgrid_circuit )
% JACOBLOADBM Loads BM from file for a provided tspan
%
%       variable = 'BMx'
%       variable = 'TMx'
%

T = Time(fMechResult.time, 's');
T = T.(unit);

[d1, t1] = min(abs(T - tspan(1)));
[d2, t2] = min(abs(T - tspan(end)));

if any([d1,d2] > mean(diff(T)))
    warning('Difference of time points too large: %g, %g\n\tmax: %g.', ...
        d1, d2, mean(diff(T)));
end

BM_displacement = 'BMx';
OHC_stereocilia = 'TMx';

BM_displ = fMechResult.(BM_displacement)(t1:t2,:);
OHC_cilia = fMechResult.(OHC_stereocilia)(t1:t2,:);
t = T(t1:t2,:);

nt = numel(t);
if nt == 1
    F_IHC = griddedInterpolant(xgrid_mech, BM_displ);
    F_OHC = griddedInterpolant(xgrid_mech, OHC_cilia);
    BM_displ_circuit = F_IHC(xgrid_circuit);
    OHC_cilia_circuit = F_OHC(xgrid_circuit);
else
    F_IHC = griddedInterpolant({t,xgrid_mech}, BM_displ);
    F_OHC = griddedInterpolant({t,xgrid_mech}, OHC_cilia);
    BM_displ_circuit = F_IHC({t,xgrid_circuit});
    OHC_cilia_circuit = F_OHC({t,xgrid_circuit});
end

vars = {t, BM_displ_circuit, OHC_cilia_circuit};

end
