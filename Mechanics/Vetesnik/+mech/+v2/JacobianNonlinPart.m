function [ J ] = JacobianNonlinPart(y, N, BMy, TMa, YCutDer)
arguments
    y         double
    N   (1,1) double
    BMy (:,1) double
    TMa (:,1) double
    YCutDer (1,1) function_handle   
end

dFdY = YCutDer(y(2*N+1:3*N));

O = zeros(N);
    
J23 = -diag(BMy.*dFdY);
J43 = diag(TMa.*BMy.*dFdY);

J = [ ...
      O   O   O   O; ...
      O   O J23   O; ...
      O   O   O   O; ...
      O   O J43   O];

end