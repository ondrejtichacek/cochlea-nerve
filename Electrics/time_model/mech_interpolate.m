function [BM_displ, OHC_cilia] = mech_interpolate(t, y, xpos, size_info, mechopt)
%MECH_INTERPOLATE 
arguments
    t double
    y double
    xpos double
    size_info (1,1) struct
    mechopt (1,1) mechOpt
end

m = size_info.mech_BMx.start;
n = size_info.mech_BMx.end;
BM_displ = mechopt.NMkonst_BM * y(m:n,:);

m = size_info.mech_TMx.start;
n = size_info.mech_TMx.end;
OHC_cilia = mechopt.NMkonst_OHC_cilia * y(m:n,:);

xgrid_mech = mechopt.xgrid;
xgrid_circuit = xpos;

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

BM_displ = BM_displ_circuit;
OHC_cilia = OHC_cilia_circuit;

end

