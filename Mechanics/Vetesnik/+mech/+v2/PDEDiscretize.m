function [x, BMx, BMv, BMy, BMin, GG, TMy, TMw, TMa, damp, M] = PDEDiscretize(mechopt, args)
% PDEDISCRETIZE
%
% Here we connect the physics description of the system to a PDE that we
% solve numerically
% 

% [1  0  0 ]  [X]      [ 0       1        0       0      ][X]   [ 0              ]   [0             ]
% [0  GG 0 ] d[V]/dt = [-BMx    -BMv      0       0      ][V] + [ BMin*inp(t)    ] + [ -BMy*S(Y)    ]
% [0  0  1 ]  [Y]      [ 0       0        0       1      ][Y]   [ 0              ]   [0             ]
% [0  0  GG]  [W]      [ TMa*BMx TMa*BMv -GG*TMy -GG*TMw ][W]   [-TMa*BMin*inp(t)]   [ -TMa*BMy*S(Y)]

arguments
    mechopt (1,1) mechOpt
    args.xgrid (1, :) = mech.gaussgrid(mechopt.Numstacks);
    args.plotflag (1,1) logical = false
    args.hfig (1,:) = [];
    args.damping_factor (1,1) double = 1.0 % for over_damping (values other than 1 for special use only)
end

if args.plotflag == true
    if isempty(args.hfig)
        args.hfig = struct('undamp', figure);
    end
end

% args.xgrid = linspace(0,1, mechopt.Numstacks+1);
% args.xgrid = args.xgrid(2:end);

x = args.xgrid;

% mechanical parameters of the TM
[gamma, omega2] = mech.v2.TM_params(x, ...
        mechopt.TM_tuning_shift, ...
        mechopt.BM_tuning, ...
        mechopt.TM_damping_coef, ...
        'plotflag', args.plotflag);

% [gamma, omega2] = mech.v1.TM_params(x, ...
%         'plotflag', args.plotflag);

% [gamma, omega2] = mech.v1.TM_params(x*0.9);

TMy = omega2;
TMw = gamma;
TMa = ones(mechopt.Numstacks,1);

if isfield(mechopt.user_settings, 'experiments')
    experiments = mechopt.user_settings.experiments;
else
    experiments = struct( ...
        'fixed_BM', false);
end

% mechanical parameters of the BM
[stiff, damp, shear, undamp, Gs, GG, G, M] = mech.v2.BM_params( ...
        gamma, x, ...
        mechopt.gain, ...
        mechopt.lambda, ...
        mechopt.shearing_coefficient, ...
        mechopt.bmdata_params, ...
        'experiments', experiments, ...
        'plotflag', args.plotflag);

damp = args.damping_factor * damp;

% TMy = TMy .* diag(M) * 2;
% TMw = TMw .* diag(M);

% figure
% tiledlayout flow
% nexttile
% hold on
% plot(TMy)
% plot(stiff)
% set(gca, 'YScale', 'log')
% nexttile
% hold on
% plot(TMw)
% plot(diag(damp))

BMx = stiff;
BMv = damp + shear;
BMy = undamp;

BMin = Gs;

if do_plot('undamp')
    figure(args.hfig.undamp);
    plot(x, undamp, 'k');
    xlabel('fractional distance from stapes')
    ylabel('BM undamping');
end

function tf = do_plot(figname)
    tf = (args.plotflag == true && isfield(args.hfig, figname));
end

end
