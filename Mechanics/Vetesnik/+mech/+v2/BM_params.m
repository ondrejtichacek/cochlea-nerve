function [ stiff, damp, shear, undamp, Gs, GG, G, M ] = BM_params( bigamma, x, gain, lambda, shearing_coefficient, bmdata_params, args )
%BM_PARAMS
% v2
% parameters of basilar membrane
arguments
    bigamma
    x      (:,1) double
    gain   (1,1) double  = 0.87
    lambda = 1.157
    shearing_coefficient (1,1) double = 95
    bmdata_params cell = {}
    args.experiments = []
    args.plotflag (1,1) logical = false
end


% OUTPUT
% Gs ... represents the ratio between the BM displacement and
%        stapes displacement as the direct effect of stapes motion
%        on the BM and the BM-to-BM reaction at zero BM stiffness
%        and viscosity (BMinput decays rapidly from base to apex
%        starting from about 250)
%
% GG ... = G+M is the BM integro-differential motion equation kernel
%        where M ... mass
%              G ... Green's function matrix G

% [GG, ~, stiff, dampn, damp0, shear, Gs, sh, G, M] = mech.v1.bmdata(x, ...
%     bmdata_params{:}, ...
%     'verflag', 'v2', 'plotflag', args.plotflag);

[GG, ~, stiff, dampn, damp0, shear, Gs, sh, G, M] = mech.v1.bmdata(x, ...
    bmdata_params{:}, ...
    'verflag', 'v3', 'plotflag', args.plotflag);

% [GG, ~, stiff, dampn, damp0, shear, Gs, sh, G, M] = mech.v1.bmdata(x, ...
%     bmdata_params{:}, ...
%     'verflag', 'v4', 'plotflag', args.plotflag);

figure
hold on
set(gca, 'Yscale', 'log');
plot(x, sqrt(stiff./diag(M))/2/pi)
[~, cf] = characteristic_frequency_model('Greenwood');
plot(x, cf(x));
[~, cf] = characteristic_frequency_model('Li');
plot(x, cf(x));


% dampn = damp0;

%%

% shearing_coefficient = 95; % Vetesnik
% shearing_coefficient = 25; % Nobili
shear = shearing_coefficient * shear;
sh = shearing_coefficient * sh;

damp0_orig = damp0;

% damp0 = damp0 - damp0(1) + dampn(1) + sh;

% undamp = gain*damp0.*bigamma.*lam;
undamp = mech.undamping(bigamma, damp0, gain, lambda, x);

% g1 = 1.44;
% g2 = 5.88;
% g3 = 0.15;
% L = 1;
% u = g1 * (1+g2*tanh(g3*(1-x/L)));
% undamp = mech.undamping(u, 1, gain, lambda, x);
% A = -200;
% undamp = A * damp0 .* u * gain;

if args.plotflag == true
    figure
    hold on
    
    plot(x, damp0_orig, 'k--', 'DisplayName','damp0_{orig}');
    plot(x, damp0, 'k', 'DisplayName','damp0');
    plot(x, undamp./bigamma, 'k:', 'DisplayName','\lambda\cdotdamp0');
    
    plot(x, dampn, 'b', 'DisplayName','dampn');
    plot(x, dampn + sh, 'b--', 'DisplayName','dampn+sh');
    
    plot(x, sh, 'r--', 'DisplayName','sh');
    plot(x, diag(shear, 0), 'r-', 'DisplayName','Sh_{0}');
    plot(x(2:end), diag(shear, -1), 'r:', 'DisplayName','Sh_{-1}');
    plot(x(1:end-1), diag(shear, 1), 'r:', 'DisplayName','Sh_{+1}');
    
    plot(x, dampn + diag(shear), 'm', 'DisplayName','BMv_{0}');

    title('Damping and undamping')
    ylabel('(un)damping')
    xlabel('fractional distance from stapes')
    legend
end

damp = sparse(diag(dampn));
    
end