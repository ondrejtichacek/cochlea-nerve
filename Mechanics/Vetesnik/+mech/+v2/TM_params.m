function [ gamma, omega2 ] = TM_params( x, x0, BM_tuning, damping_coef, args )
%TM_PARAMS
% v2
% Tuning parameters of the tectotiral membrane
arguments
    x (:,1) double = linspace(0,1,300)
    x0 (:,1) double = 0.05
    BM_tuning (1,1) function_handle = @(x) 165.4 * (10.^(2.1*(1-x)) - 0.88);
    damping_coef (1,1) double = 21.50;
    args.plotflag (1,1) logical = false
end

% resonance_frequency = BM_tuning(x + x0) .* fac;
% see eg. https://www.pnas.org/doi/10.1073/pnas.93.16.8727

fac = 0.75; % 0.75 correspond to half-octave shift,
% But this shift creates a large CF offset  for low frequencies

fac = linf(0.75, 0.05, x);
% fac = 0.75 + boltzmann(x, -0.7, 0.8, 0.02);
% fac = 0.75 + boltzmann(x, 0.5, 0.8, 0.02);
resonance_frequency = BM_tuning(x + x0) .* fac;

% resonance_frequency = BM_tuning(x + x0);
resonance_frequency = max(resonance_frequency, 5.0);

assert(all(resonance_frequency >= 5.0));

omega = 2*pi * resonance_frequency;

% omega2 is equivalent to K/M along the TM
%     K ... stiffness coef
%     M ... mass
omega2 = omega.*omega;

% zeta = damping_coef;
% zeta = damping_coef / 33;
zeta = damping_coef ./ sqrt(resonance_frequency);

% zeta = zeta * 10;

% gamma is equivalent to H/M along the TM
%     H ... viscosity coef
%     M ... mass
gamma = 2 * zeta .* omega;

% omega_BM = BM_tuning(x);
% gamma = 2 * zeta .* omega_BM;

N = numel(x);
assert(all(size(gamma) == [N,1]));
assert(all(size(omega2) == [N,1]));

if args.plotflag == true
    figure
    hold on
    plot(x, BM_tuning(x), 'b')
    plot(x, resonance_frequency, 'r')
    xlabel('fractional distance from stapes')
    ylabel('Frequency (Hz)')
    legend('Basilar Membrane', 'Tectorial Membrane')
    title('Resonance frequency');
    
    figure
    plot(x, omega2, 'k')
    xlabel('fractional distance from stapes')
    ylabel('\omega_{TM}^2 [rad^2/s^2]')
    
    figure
    plot(x, zeta, 'k')
    xlabel('fractional distance from stapes')
    ylabel('\zeta_{TM}')
    
    figure
    hold on
    plot(x, omega, 'b', 'DisplayName', '\omega_{TM}' )
    plot(x, gamma, 'r', 'DisplayName', '\gamma_{TM}')
    xlabel('fractional distance from stapes')
    ylabel('[rad/s]')
    legend()
end

end