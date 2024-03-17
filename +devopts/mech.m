function [midopt, mechopt] = mech(GlobalSamplingFrequency, args)
arguments
    GlobalSamplingFrequency (1,1) Frequency = Frequency(200, 'kHz');
    
    args.Numstacks (1,1) double = 600;
    
    args.amplifier (1,:) char = 'electric'
    % args.amplifier (1,:) char = 'mechanic'
    args.gain_factor (1,1) double = 1;
    
    args.opt_ver (1,:) char = 'latest'
    args.plotopt (1,:) plotOpt = plotOpt('MNA');
    
    args.oe_identifier (1,:) char = 'Meddis'
    args.me_identifier (1,:) char = 'PBLL'
    
    args.noise_damage (1,1) logical = false;
    args.noise_damage_f0 (1,1) Frequency = Frequency(1200, 'Hz')

    args.plotflag = false;
end

plotopt = args.plotopt;

%% Outer Ear options

% The model can run with or without OE.
% At the moment, OE is ignored when ME is set to 'none'

switch args.oe_identifier
    case 'Meddis'

        externalResonanceFilters = [ ... gain order lp hp
          -11.99927987142523,  2,   206.2119974669669,   358.5102170379319;
          -11.99998556203806,  2,   288.0827584566169,   444.8320921474029;
          -6.660289825314575,  2,    604.601203475762,   1656.284925736281;
          -6.009994624924989,  2,   365.6631053098679,   1177.111851175218;
          -4.392649142756406,  2,   1600.417554322875,   3593.838940335973;
          -3.819930936639489,  2,   1634.013499558196,   3583.542619484177;
         -0.9495775823395274,  2,   2415.056217002997,   5319.282950031623;
           8.077026047557176,  2,   2435.063884295232,   4541.144792545283;
           3.016053335238496,  2,   2809.307899606387,    7618.45593491728;
           6.209711798089334,  2,   1398.641611236509,   4696.586038747771;
          -9.118556974337517,  2,   6461.689628623291,   9721.286323027263;
        ];
        
        outer_ear_params = struct( ...
                'identifier', 'Meddis', ...
                'args', struct( ...
                    'externalResonanceFilters', externalResonanceFilters));
    case 'none'
        outer_ear_params = struct( ...
                'identifier', 'none');
    otherwise
        error('Not set up')
end
        
%% Middle Ear options

midopt = midOpt.simple_init('dev');

% The model can run with or without ME.
% In case of no middle ear, args.me_identifier = 'none'
% In that case, the mock middle ear has a flat frequency response,
% comparable to the results with ME at around 1 kHz. Consequently, 
% vd_OW is chosen such that with respect to AMo and L the stapes
% acceleration for a 1 kHz pure tone is approx. equal:
%  - no ME:   395
%  - ME:      390
%  - ME + OE: 351
% For other frequencies, this will differ (as intended)

switch args.me_identifier
    case 'PBLL'
        me_params = struct(...
            'identifier', 'PBLL', ...
            'vd_OW', 3.56);

    case 'none'
        me_params = struct(...
            'identifier', 'none', ...
            'AMo', 1e-5, ...
            'L', 0);
    otherwise
        error('Not set up')
end

%% Mech options

if strcmp(args.amplifier, 'electric')
    
    % GAIN = 30;
    GAIN = 1; % for a steeper MET transduction function
    % GAIN = 0;
    % GAIN = 0.9;
    
    % theta should be scaled in ~ 10-50 mV => 1e-3
    % theta = 1e-3 * 25; % mV


    switch args.opt_ver
        case {'latest', 'm6'}

            cxx = [1, 0.9, 0.8, 0.7, 0.6, 0.2, 0.1, 0];
            cx = cxx;

            lambda = 1.0 * [
                1.5
                1.5
                1.4 %1.4 %1.3 %1.2 %1.1
                1.4
                1.5
                1.7 %1.6 %1.5 %1.1
                1.9 %1.7 %1.6 %1.5 %1.0
                2.0 %1.8 %1.7 %1.5 %1.0
            ]';

%             cxx = [1, 0.9, 0.8, 0.2, 0.1, 0];
%             cx = cxx;
% 
%             lambda = 1.0 * [
%                 1.5
%                 1.5
%                 1.2 %1.1
%                 1.5 %1.1
%                 1.6 %1.5 %1.0
%                 1.7 %1.5 %1.0
%             ]';

            theta = ones(size(lambda));
        case {'m5'}

            cxx = [1, 0.81, 0.61, 0.45, 0.33, 0.21, 0.09, 0];
            cx = cxx;

            lambda = 1.0 * [
                1.1
                1.1
                1.2
                1.3
                1.5
                1.4
                1.3
                0.9
            ]';

            theta = ones(size(lambda));

        case {'m4'}

            cxx = [1, 0.81, 0.61, 0.45, 0.33, 0.21, 0.09, 0];
            % cx = [1, cxx, 0];
            cx = cxx;

            lambda = 1.0 * [
            ...1.1976997423422096
            1.1
            ...1.2340890901729220
            1.1
            0.8613146206062111
            1.3436661539907031
            1.4451903498387002
            1.3190878589270696
            1.2375900930971571
            0.8462296865113116
            ]';

            theta = ones(size(lambda));

        case {'m3'}

            cxx = [1, 0.81, 0.61, 0.45, 0.33, 0.21, 0.09, 0];
            % cx = [1, cxx, 0];
            cx = cxx;

            lambda = 1.0 * [
            0.9508148090066133
            1.1736556648067573
            1.0120781918498101
            1.3832596454440320
            1.4235736378927197
            1.2087196700894809
            1.2867387552823857
            0.9409488180613454
            ]';

            theta = ones(size(lambda));

        case {'m2'}
            
            cxx = [0.89 0.81 0.69 0.61 ...
                   ... 0.57 0.49 0.45 0.41 0.37 0.33 0.29 0.25 0.21 0.17 0.13 0.09 
                   0.05];
            cx = [1, cxx, 0];

            lambda = 1.0 * [
            0.80      ... 1
            ...0.90      ... 0.93
            0.80      ... 0.89
            0.80      ... 0.81
            1.00      ... 0.69
            1.00      ... 0.61
            ...1.00      ... 0.57
            ...1.00      ... 0.49
            ...1.00      ... 0.45
            ...1.00      ... 0.41
            ...1.00      ... 0.37
            ...1.00      ... 0.33
            ...1.00      ... 0.29
            ...1.00      ... 0.25
            ...1.00      ... 0.21
            ...1.00      ... 0.17
            ...1.00      ... 0.13
            ...1.00      ... 0.09
            1.70      ... 0.05
            1.70      ... 0
            ];

            % lambda = lambda -0.15 + 0.3 * (1-cx);

            theta = ones(size(lambda));

        otherwise
            error('Unknown version %s', args.opt_ver)
    end
end
   

if strcmp(args.amplifier, 'mechanic')
    
    GAIN = 1;
    
    % FREQUENCY = [199, 328, 495, 712, 996, 1364, 1843, 2466, 3276, 4330, 5701, 7485];
    % cxx = [0.8172, 0.7296, 0.6575, 0.5938, 0.5350, 0.4799, 0.4272, 0.3761, 0.3264, 0.2775, 0.2293, 0.1816];
    % cx = [1, cxx, 0];
    % cxx = [0.93 0.81 0.69 0.57 0.45 0.41 0.37 0.33 0.29 0.25 0.21 0.17 0.13 0.09 0.05];
    
    cxx = [0.69, 0.57, 0.45, 0.41, 0.37, 0.33, 0.29, 0.25, 0.21];
    cx = sort([1, 0.93, 0.81, ...
            cxx, ...
          0.17, 0.13, 0.09, 0.05, 0], 'desc');
    
    % mt = struct( ...
    %     'base', 2.4531, ...
    %     'alpha', -6.36, ...
    %     'Fo', 21100);
    % 
    % cf = @(x) mt.Fo * mt.base.^(mt.alpha*(x));
    % cx = @(c) log(c/mt.Fo)/log(mt.base^mt.alpha);

switch args.opt_ver
    case 'latest'
        lambda = [ % opt_14d
            0.6
            0.6
            0.6
        0.8339745460527473
        1.0545264935514496
        1.1044090703494349
        1.0607643125062187
        1.2171054558286960
        1.0530045036156734
        1.1350302231501419
        1.2621714506079365
        0.7087437701453470
            0.95
            0.95
            0.9
            1.0
            1.0
        ];
    otherwise
        error('Unknown version %s', args.opt_ver)
end
switch args.opt_ver
    case 'latest'
        theta = 1e-3 * [ % opt_14d
            1500.0
            1500.0
            1500.0
        1999.3373835245149621
        56.1885752530849558
        76.3014204687202664
        102.7857880984287817
        31.4624726019959304
        705.9802714701031618
        9.4487591209666704
        74.2861408671076902
        1194.2599641840943150
            1000.0
            1000.0
            1000.0
            1000.0
            1000.0
        ];
    otherwise
        error('Unknown version %s', args.opt_ver)
end

end

%%

if args.plotflag == true
    hfig = plotopt.figure();
    tiledlayout(4,1)
else
    hfig = [];
end

if args.plotflag == true
    nexttile
end

lambda_gain = interpolateAlongBM( ...
         lambda, ...
    'x', cx, ...
        'plotflag', args.plotflag, ...
        'hfig', hfig);
if args.plotflag == true
    ylabel('\lambda gain magnitude')
end

%%
if args.noise_damage

    [cxfun, cf] = characteristic_frequency_model('Li');
    % [cxfun, cf] = characteristic_frequency_model('Greenwood');

    % noise_band = [800, 1600];
    noise_band = [args.noise_damage_f0.Hz * 2/3,  args.noise_damage_f0.Hz * 4/3];
    noise_band = sort(noise_band, 'descend');
    
    % v = 3; % only set up for 800-1600 noise band
    v = 1;
    v = 1.1;
    if v == 1
        noise_band_cx = cxfun(noise_band);
    
        % D = [-0.05, 0.02];
        % D = [-0.1, 0.15];
        % D = [-0.025, 0.1];% v1_1
        D = [-0.05, 0.1];

        % d = 0.001 * [1, -1];
        % d = [-0.05, 0.1];
        % d = [0.025, 0.05];% v1_1
        d = [0.00, 0.05];
    
        noise_band_cx_outer = noise_band_cx + D;
        noise_band_cx_inner = noise_band_cx + d;

        lambda_outer = lambda_gain(noise_band_cx_outer); % orig gain
        lambda_inner = [0, 0];

    elseif v == 1.1
        noise_band_cx = cxfun(noise_band);
    
        D = [-0.05, 0.1];

        d = [0.00, 0.05];
    
        noise_band_cx_outer = noise_band_cx + D;
        noise_band_cx_inner = noise_band_cx + d;

        lambda_outer = lambda_gain(noise_band_cx_outer); % orig gain
        lambda_inner = [0, 0];

    elseif v == 2
        % manually selecting the area - performing simulations of high SPL
        % tones of 800 and 1600 Hz

        noise_band_cx_outer = [0.35, 0.8];
        noise_band_cx_inner = [0.45, 0.75];
    
        lambda_outer = lambda_gain(noise_band_cx_outer); % orig gain
        lambda_inner = [0, 0];

    elseif v == 3
        % manually selecting the area - performing simulations of high SPL
        % tones of 800 and 1600 Hz

        noise_band_cx_outer = [0.4, 0.8];
        noise_band_cx_inner = [0.5, 0.75];
    
        lambda_outer = lambda_gain(noise_band_cx_outer); % orig gain
        lambda_inner = [0, 0];
    end
    
    ind = cx >= noise_band_cx_outer(1) & cx <= noise_band_cx_outer(2);
    
    lambda_new = lambda;
    cx_new = cx;

    lambda_new(ind) = [];
    cx_new(ind) = [];
    
    lambda_new = [lambda_new, lambda_outer, lambda_inner];
    cx_new = [cx_new, noise_band_cx_outer, noise_band_cx_inner];
    
    lambda_gain_orig = lambda_gain;

    lambda_gain = interpolateAlongBM( ...
             lambda_new, ...
        'x', cx_new, ...
            'plotflag', args.plotflag, ...
            'hfig', hfig);

    % Extra plot
    hfig = figure;
    xxx = linspace(0,1,256);
    plot(xxx, lambda_gain(xxx) ./ lambda_gain_orig(xxx), ...
        Color=linecolors(5), ...
        DisplayName='Amplifier strength')

    yl = ylim();
    
    ylim(yl - [0.1, 0]);

    h = ShadePlotForEmpahsis(noise_band_cx_inner, 'k', 0.1);
    h.DisplayName = 'Noise band';
    legend('Location', 'southeast')

    ylabel('Rel. amplif.')
    xlabel('Normalized distance form stapes')

    ylim(yl - [0.1, 0]);
end
%%

% nonlinearity_sensitivity of value \theta means that the sigmoidal
% function is:
%    ~linear    [-0.2 to  0.2]
%    saturated  [-inf to -0.5]
%    saturated  [ 1.5 to  inf]
% the function is not symmetric.

if args.plotflag == true
    nexttile
end

nonlinearity_sensitivity = interpolateAlongBM( ...
         theta, ...
    'x', cx, ...
        'plotflag', args.plotflag, ...
        'hfig', hfig);
if args.plotflag == true
    ylabel('\theta gain saturation')
end

x = mech.gaussgrid(args.Numstacks);
% x = linspace(0, 1, args.Numstacks+1);
% x = x(2:end);

% nonlinearity_sensitivity = @(x) 1;

amplifier = amplifier_params( ...
    nonlinearity_sensitivity(x), ...
    'hfig', hfig);

if args.plotflag == true
    nexttile([2 1])
    
    for i = 1:numel(cx)
        [~, xx(i)] = find(x >= cx(i), 1);
    end
    xx = sort(xx, 'desc');
    c = rwb(numel(xx));
    M = max(nonlinearity_sensitivity(x));
    y = linspace(-M, M, 200);
    
    hold on
    for i = 1:numel(xx)
        Y = nonlin(y, amplifier.y1(xx(i)), amplifier.y2(xx(i)), amplifier.c1, amplifier.c2, amplifier.b(xx(i)), amplifier.q(xx(i)));
        %Y = lambda_gain(x(xx(i))) * Y;
        plot(y*1e3, Y(:)./y(:)*100, 'Color', c(i,:), 'DisplayName', sprintf('x = %g', xx(i)/args.Numstacks));
    end
    legend();
    xlabel('rel. Voltage (mV)');
    ylabel('rel. amplification (%)');
end


if args.plotflag == true
    hfig = plotopt.figure();
    tiledlayout(2,1)
    nexttile


    amplifier_params( ...
        nonlinearity_sensitivity(x), ...
        'X_plot_scale', 1000, ...
        'Y_plot_scale', 1000, ...
        'plotflag', args.plotflag, ...
        'hfig', hfig);


    xlabel('rel. Voltage (mV)');
    ylabel('Response (a.u.)');
end
    
if args.plotflag == true
    nexttile
    
    for i = 1:numel(cx)
        [~, xx(i)] = find(x >= cx(i), 1);
    end
    xx = sort(xx, 'desc');
    c = rwb(numel(xx));
    M = max(nonlinearity_sensitivity(x));
    y = linspace(-M, M, 200);
    
    hold on
    for i = 1:numel(xx)
        Y = nonlin(y, amplifier.y1(xx(i)), amplifier.y2(xx(i)), amplifier.c1, amplifier.c2, amplifier.b(xx(i)), amplifier.q(xx(i)));
        %Y = lambda_gain(x(xx(i))) * Y;
        plot(y*1e3, Y*1e3, 'Color', c(i,:));
    end
    
    xlabel('rel. Voltage (mV)');
    ylabel('Response (a.u.)');
end
    
%%

% TM_damping_coef ... the higher value, the more TM is damped
% shearing_coeficient ... relative viscosity of OC compared to water

mechopt = mechOpt( ...
    'middle_ear_params', me_params, ...
    'outer_ear_params', outer_ear_params, ...
    ...
    'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
    'Numstacks', args.Numstacks, ...
    ...
    'specifier', ...
        ...'dev', ...
        'v2', ...
        ...'v1', ...
        ...'v0', ...
    ...
    ...'Solver', @odeEuler, ...
    ...'Solver', @odeDormandPrince, ...
    ...'Solver', @BDF2Solver, ...
    'Solver', @ode15s, ...
    ...'Solver', @ode45, ...
    'NumDiv', 'auto', ...
    ...'save_method', 'c_posix', ...
    ...'save_method', 'matlab_fwrite', ...
    'save_method', 'matlab_matfile', ...
    ...
    ...'amplifier', 'none', ...
    ...'amplifier', 'mechanic', ...
    ...'amplifier', 'electric', ...
    'amplifier', args.amplifier, ...
    ...
    ...'lambda', 0, ...
    ...'lambda', 1.157, ...
    ...'lambda', '+mech/lambdas6a.mat', ...
    ...'lambda', '+mech/lambdas6.mat', ...
    'lambda', lambda_gain, ...
    ...
    'nonlinParams', amplifier, ... % Dev
    ...'nonlinParams', amplifier_params(0.2), ... % Dev
    ...'nonlinParams', amplifier_params(1), ... % Nobili
    ...'nonlinParams', amplifier_params(0.1), ... % Vetesnik
    ...
    'identifier', ...
        'Vetesnik' ...
    ...    'Verhulst' ...
    );

if mechopt.amplifier == "mechanic"

    mechopt.approximation = 'nonlinear';

    mechopt.TM_tuning_shift = 0.05;
    mechopt.TM_damping_coef = 21.5;
    mechopt.shearing_coefficient = 95;
    
    mechopt.NMkonst_BM = 20;
    mechopt.NMkonst_OHC_cilia = 20;

    % mechopt.gain = 0.95;
    mechopt.gain = GAIN * args.gain_factor;

    % mechopt.integration = 'standalone';
    mechopt.integration = 'electric';

elseif mechopt.amplifier == "electric"

    mechopt.approximation = 'linear';

    % KK = 1;
    % KK = 3;
    KK = 6;
    switch args.opt_ver
        case 'm4'
            KK = 6;
    end

    mechopt.TM_tuning_shift = 0.25;
    % mechopt.TM_tuning_shift = 0.2336023268388646;
    % mechopt.TM_tuning_shift = 0.1806523586396664;
    % mechopt.TM_tuning_shift = 0.1806523586396664; % v0
    % mechopt.TM_tuning_shift = 0.1833702249636318; % opt4a
    % mechopt.TM_tuning_shift = 0.2576674459168363; % opt5a
    % mechopt.TM_tuning_shift = 0.2624957506332676; % opt5b
%     mechopt.TM_tuning_shift = 0.2532742098171023; % opt5c
    
    % mechopt.TM_tuning_shift = 0.01;
    % mechopt.TM_tuning_shift = 0.05;
    % mechopt.TM_tuning_shift = 0.00;
    mechopt.TM_tuning_shift = 0.05;

    % mechopt.TM_damping_coef = 21.5;
    % mechopt.TM_damping_coef = 50;
    % mechopt.TM_damping_coef = 54.8503040069368311;
    % mechopt.TM_damping_coef = 45.7681247909224922;
    % mechopt.TM_damping_coef = 30; % v0
    % mechopt.TM_damping_coef = 30.5546955820619175; % opt4a
    % mechopt.TM_damping_coef = 22.7782550569935438; % opt5a
    %mechopt.TM_damping_coef = 23.6501028155006026; % opt5b
    mechopt.TM_damping_coef = 26.0801183314589764; % opt5c
    
    mechopt.shearing_coefficient = 95;
    mechopt.shearing_coefficient = 195;
%     mechopt.shearing_coefficient = 20.9;
    
    mechopt.NMkonst_BM = 20;
    mechopt.NMkonst_OHC_cilia = 20/KK;

    mechopt.prestin_voltage_sensitivity = - KK*20;
    % mechopt.prestin_voltage_sensitivity = KK*3.5;

    switch args.opt_ver
        case 'm4'
            mechopt.prestin_voltage_sensitivity = -100;
            % mechopt.prestin_voltage_sensitivity = -KK * 55;
            % mechopt.prestin_voltage_sensitivity = -KK * 48.3776978010050058
            % mechopt.prestin_voltage_sensitivity = -KK * 49.5719783206995359;
            % mechopt.prestin_voltage_sensitivity = -KK * 50; % v0
            % mechopt.prestin_voltage_sensitivity = -200.9115621775761156;
            % mechopt.prestin_voltage_sensitivity = -277.3378572899809456; % opt5a
            % mechopt.prestin_voltage_sensitivity = -271.4961702668460930; % opt5b
            % mechopt.prestin_voltage_sensitivity = -269.7843481351479227; % opt5c
    end


    mechopt.bmdata_params = { ...
        ...'stiff_exp', -3.6, ...
        ...'stiff_0', 2000, ...
        ...'stiff_exp', -3.1833, ...
        ...'stiff_0', 4529, ...
        ...'stiff_exp', [-3.1, -1.0, -5], ...
        ...'stiff_0', [4500, -5, -2000], ...
        'stiff_exp', [-2.8, -1.0], ...
        'stiff_0', [3500, -55], ...
        ...'stiff_exp', [-4.2], ...
        ...'stiff_0', [2000], ...
        ...'stiff_exp', [-2.566], ...
        ...'stiff_0', [2000], ...
        };

    mechopt.gain = 1 * args.gain_factor;
    % mechopt.gain = GAIN * args.gain_factor;

    % mechopt.integration = 'standalone';
    mechopt.integration = 'electric';

else
    error('Something is wrong');
end

% to be adjusted in the future
if mechopt.gain == 0
    mechopt.char_frequency_model_best_fit = 'Greenwood';
else
    mechopt.char_frequency_model_best_fit = 'Li';
end

% mechopt = mechOpt( ...
%     'samplingFrequency', GlobalSamplingFrequency, ... [Hz]
%     'Numstacks', 1000, ...
%     'tf', tf, ...                [s]
%     'irregularities', [1], ...
%     'normalizeRMS', [1], ...
%     'subjectNo', 1, ...
%     'MechModel', ...
%     ...'Vetesnik' ...
%     'Verhulst' ...
%     );

end