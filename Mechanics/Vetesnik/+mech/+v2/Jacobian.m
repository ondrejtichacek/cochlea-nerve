function [ J ] = Jacobian( N, GG, BMx, BMv, BMy, TMy, TMw, TMa, approximation, amplifier)
arguments
    N   (1,1) double
    GG  (:,:) double
    BMx (:,1) double
    BMv (:,:) double
    BMy (:,1) double
    TMy (:,1) double
    TMw (:,1) double
    TMa (:,1) double
    approximation char = 'nonlinear'
    amplifier char = 'mechanic'
end

I = eye(N);
O = zeros(N);

J21 = -diag(BMx);
J22 = -BMv;

J41 = diag(TMa.*BMx);
J42 = diag(TMa)*BMv;
J43 = -GG*diag(TMy);
J44 = -GG*diag(TMw);
% J43 = -diag(TMy);
% J44 = -diag(TMw);

J = [ ...
      O   I   O   O; ...
    J21 J22   O   O; ...
      O   O   O   I; ...
    J41 J42 J43 J44];

% J = [ ...
%       O   I   O   O; ...
%     J21 J22   O   O; ...
%       O   O   O   I; ...
%       O   O J21 J22];


if strcmp(approximation, 'linear') && strcmp(amplifier, 'mechanic')
    
    J23 = -diag(BMy);
    J43 = diag(TMa.*BMy);

    J = J + [ ...
              O   O   O   O; ...
              O   O J23   O; ...
              O   O   O   O; ...
              O   O J43   O];
end

end

