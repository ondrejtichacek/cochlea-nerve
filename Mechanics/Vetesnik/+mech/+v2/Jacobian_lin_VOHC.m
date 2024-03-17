function [ J2, J4 ] = Jacobian_lin_VOHC(BMy, TMa, const_OHC_force)
arguments
    BMy (:,1) double
    TMa (:,1) double
    const_OHC_force (1,1) double
end

dFdY = const_OHC_force;
% dFdY = - const_OHC_force;
    
J2 = -BMy.*dFdY;
J4 = TMa.*BMy.*dFdY;

end