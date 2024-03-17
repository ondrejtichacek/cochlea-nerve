function [ J2, J4 ] = JacobianNonlinPart_VOHC(vohc, N, BMy, TMa, YCutDer, const_OHC_force)
arguments
    vohc      double
    N   (1,1) double
    BMy (:,1) double
    TMa (:,1) double
    YCutDer (1,1) function_handle
    const_OHC_force (1,1) double
end

dFdY = const_OHC_force * YCutDer(const_OHC_force * vohc);
% dFdY = - const_OHC_force * YCutDer(const_OHC_force * vohc);
    
J2 = -BMy.*dFdY;
J4 = TMa.*BMy.*dFdY;

end