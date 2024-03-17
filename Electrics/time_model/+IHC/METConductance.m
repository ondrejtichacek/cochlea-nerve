function [ ...
    MET_conductance, MET_conductance_der, ...
    MET_conductance_pre, MET_conductance_der_pre, ...
    MET_popen_BM_mag_pre, G_max_pre] = METConductance( xgrid, simulation_units, args )
%METCONDUCTANCE
%IHC

arguments
    xgrid (:,1) double = linspace(0, 1, 300)'
    simulation_units (1,1) struct = struct('voltage', 'V')
    args.version = 'Jia_2007'
    args.plotflag = false
    args.save_figures = false
    args.handles = struct()
    args.optim_params = []

    args.noise_damage (1,1) logical = false
    args.noise_damage_f0 (1,1) Frequency = Frequency(1200, 'Hz')
    args.noise_damage_degree (1,1) double = 0.1
end

%%

% [char_frequency, char_position] = frequency_map();
[char_position, char_frequency] = characteristic_frequency_model('Li');

default_tikz_options = { ...
    'relativeDataPath', 'img/HairCells/', ...
    'extraaxisoptions', {'legend style={font=\tiny}'}, ...
    'parseStrings', true};

%%

switch args.version

    case 'Tichacek_2022'
        
        % -----------------------------------------------------------------
        % IHC MET open probability

        % https://www.jneurosci.org/content/27/5/1006
        %
        % The approximate recording locations were 2.7–3.2 and 8.6–8.9 mm 
        % from the base of the cochlea, representing 18–13 and 1–0.6 kHz of
        % the best frequencies in vivo (Müller, 1996).
        %
        % Curve-fit parameters for the plots shown are for the basal cell 
        % and are as follows:
        %   imax = 2031.9 pA, x1 = 35.3 nm, α = 16.1 nm.
        % For the apical cell, the curve-fit parameters are as follows:
        %   imax = 2195.2 pA, x1 = 82.1 nm, α = 34.4 nm.

        x1 = (char_position(18e3) + char_position(13e3))/2;
        x2 = (char_position(1e3) + char_position(0.6e3))/2;

        % we assume linear dependence of the boltzmann parameters with
        % position perform linear regression
        % y = b0 + b1 * x
        x_jia2007 = [x1; x2];
        X = [ones(length(x_jia2007),1), x_jia2007];

        % boltzmann half-activatiobn parameter
        y_h = [35.3; 82.1];
        [h0, h1] = linear_regression(X, y_h);

        % boltzmann slope parameter
        y_s = [16.1; 34.4];
        [s0, s1] = linear_regression(X, y_s);

        MET_popen_BM_v1 = @(x,y) boltzmann(y, 1, h1*x + h0, s1*x + s0);
        MET_popen_BM_der_v1 = @(x,y) boltzmann_der(y, 1, h1*x + h0, s1*x + s0);
        
        if args.plotflag == true
            if isempty(args.handles) || numel(fieldnames(args.handles)) == 0
                hfig = figure();
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on
            
            u = linspace(-100e-9, 200e-9, 200);

            plot(u*1e9, MET_popen_BM_v1(zeros(size(u)),u*1e9), 'DisplayName', 'Jia 2007 base');
            plot(u*1e9, MET_popen_BM_v1(ones(size(u)),u*1e9), 'DisplayName', 'Jia 2007 apex');
            xlabel('displacement (nm)')
            ylabel('open probability')
        end

        %% Verhulst 2018 for comparison
        
        if args.plotflag == true
            if isempty(args.handles) || numel(fieldnames(args.handles)) == 0
                hfig = figure();
                args.handles.popen = hfig;
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on

            s0 = 40;
            u0 = 20;
            
            s1 = 16;
            u1 = 20;
            
            popen = @(u) 1 ./ (1 + exp((u0-u)/s0).*(1 + exp((u1-u)/s1)));

            u = linspace(-100e-9, 200e-9, 200);

            plot(u*1e9, popen(u*1e9), 'DisplayName', 'Verhulst 2018');
            xlabel('displacement (nm)')
            ylabel('open probability')

        end

        %% Shift
        % -----------------------------------------------------------------
        
        % Start with the same parameters as Verhulst 2018

        s0 = 40;
        u0 = 20;
        
        s1 = 16;
        u1 = 20;

        % Johnson 2011 pp 1149
        % As with OHCs, perfusing 0.02 mM Ca2+ increased the peak size of
        % the MT current and also the fraction activated at rest (Figures 
        % 8A and 8B). The mean MT current increased from
        %    0.79 ± 0.07 nA (1.3 Ca2+; n = 4)
        % to 1.72 ± 0.12 nA (0.02 Ca2+; n = 5, T = 23?C)
        % and the fraction on at rest increased from 
        %    0.045± 0.004 (1.3 Ca2+) 
        % to 0.17 ± 0.03 (0.02 Ca2+).
        % Using the latter fraction and correcting the standing current to
        % 36?C yields a resting MTconductance of 5.7 nS.
        
        % My reasoning:
        % We want to shift boltzmann to be at rest (y=0) at p=0.17.

        o = 5.82;

        u0 = u0 - o;
        u1 = u1 - o;
        
        % Now do optimization so that it fits Johnson 2015

        if isempty(args.optim_params)
            q0 = 1.0;
            q1 = 4;
            r0 = 0.47;
            r1 = 0.3;

%             q0 = 1.0;
%             q1 = -10;
%             r0 = 0.47;
%             r1 = 1.3;
        else
            q0 = args.optim_params(1);
            q1 = args.optim_params(2);
            r0 = args.optim_params(3);
            r1 = args.optim_params(4);
        end

        qq = @(x) q1*x - q0;
        rr = @(x) r1*x + r0;

%         popen_rest_target = [0.0497, 0.1667]
%         popen_rest = @(x) boltzmann3(0, 1, u0 - qq(x), u1 - qq(x), s0*rr(x), s1*rr(x));
%         popen_rest(0)
%         popen_rest(1)

        MET_popen_BM_mag = @(x,y,m) boltzmann3(y, m, u0 - qq(x), u1 - qq(x), s0*rr(x), s1*rr(x));
        MET_popen_BM_mag_pre = @(y,m) boltzmann3(y, m, u0 - qq(xgrid), u1 - qq(xgrid), s0*rr(xgrid), s1*rr(xgrid));

        MET_popen_BM_mag_der = @(x,y,m) boltzmann3_der(y, m, u0 - qq(x), u1 - qq(x), s0*rr(x), s1*rr(x));
        MET_popen_BM_mag_der_pre = @(y,m) boltzmann3_der(y, m, u0 - qq(xgrid), u1 - qq(xgrid), s0*rr(xgrid), s1*rr(xgrid));


        if args.plotflag == true
            if isempty(args.handles) || isempty(fieldnames(args.handles))
                hfig = figure();
                args.handles.popen = hfig;
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on
            
            u = linspace(-100e-9, 200e-9, 200);

            xx = zeros(size(u));
            plot(u*1e9, MET_popen_BM_mag(xx,u*1e9,ones(size(u))), 'DisplayName', 'Tichacek 2022 base');

            xx = ones(size(u));
            plot(u*1e9, MET_popen_BM_mag(xx,u*1e9,ones(size(u))), 'DisplayName', 'Tichacek 2022 apex');

            xlabel('displacement (nm)')%MET channel max conductance
            ylabel('open probability')
            legend()

        end

        % -----------------------------------------------------------------
        % IHC MET conductance
        
        % Maximal MET conductance is 29.8 ± 0.26 and 28.9 ± 0.65 nS 
        % for the apical- and basal-turn IHCs, respectively.
        
        % maximal conductance parameter
        y = [28.9; 29.8];

        % y = y * 1.6765; % to achieve resting conductance of 5.7 nS
        y = y * 1.1400; % to achieve resting conductance of 5.7 nS

        [g0, g1] = linear_regression(X, y);
        
        MET_conductance_max = @(x) g1*x + g0;

        G_max_pre = MET_conductance_max(xgrid);
        
        FAC = 1;
        
        if args.noise_damage
            [damage_fac, damage_fun] = noise_damage(args.noise_damage_f0, xgrid(:)', args.noise_damage_degree);
                
            G_max_pre = G_max_pre .* damage_fac(:);

            MET_conductance = @(x,y) MET_popen_BM_mag(x, FAC * y, (g1*x + g0) .* damage_fun(x));
            MET_conductance_der = @(x,y) FAC * MET_popen_BM_mag_der(x, FAC * y, (g1*x + g0) .* damage_fun(x));
        else
            MET_conductance = @(x,y) MET_popen_BM_mag(x, FAC * y, g1*x + g0);
            MET_conductance_der = @(x,y) FAC * MET_popen_BM_mag_der(x, FAC * y, g1*x + g0);
        end

        MET_conductance_pre = @(y) MET_popen_BM_mag_pre(FAC * y, G_max_pre);        
        MET_conductance_der_pre = @(y) FAC * MET_popen_BM_mag_der_pre(FAC * y, G_max_pre);

        if args.plotflag
            hfig = figure;
            hold on
            plot(xgrid, MET_conductance_max(xgrid), 'DisplayName', 'maximal conductance')
            plot(xgrid, MET_conductance(xgrid, zeros(size(xgrid))), 'DisplayName', 'resting conductance')
            title('IHC MET conductance')
            legend()
            xlabel('relative distance from stapes')
            ylabel('conductance (nS)')
        end

    case 'Verhulst_2018'
        
        s0 = 40e-9;
        u0 = 20e-9;
        
        s1 = 16e-9;
        u1 = 20e-9;
        
        popen = @(u) 1 ./ (1 + exp((u0-u)/s0).*(1 + exp((u1-u)/s1)));

        if args.plotflag == true
            if isempty(args.handles)
                hfig = figure();
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on
            
            u = linspace(-100e-9, 200e-9, 200);

            plot(u*1e9, popen(u), 'DisplayName', args.version);
            xlabel('displacement (nm)')
            ylabel('open probability')

        end

    case 'LopezPoveda_2006'
        
        s0 = 63.1e-9;        % Displacement sensitivity (1/m)
        u0 = 52.7e-9;        % Displacement offset (m)
        
        s1 = 12.7e-9;        % Displacement sensitivity (1/m)
        u1 = 29.4e-9;        % Displacement offset (m)
        
        popen = @(u) 1 ./ (1 + exp((u0-u)/s0).*(1 + exp((u1-u)/s1)));

        if args.plotflag == true
            if isempty(args.handles)
                hfig = figure();
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on
            
            u = linspace(-100e-9, 200e-9, 200);

            plot(u*1e9, popen(u), 'DisplayName', args.version);
            xlabel('displacement (nm)')
            ylabel('open probability')

        end

    case 'Jia_2007'
        
        % -----------------------------------------------------------------
        % IHC MET open probability

        % https://www.jneurosci.org/content/27/5/1006
        %
        % The approximate recording locations were 2.7–3.2 and 8.6–8.9 mm 
        % from the base of the cochlea, representing 18–13 and 1–0.6 kHz of
        % the best frequencies in vivo (Müller, 1996).
        %
        % Curve-fit parameters for the plots shown are for the basal cell 
        % and are as follows:
        %   imax = 2031.9 pA, x1 = 35.3 nm, α = 16.1 nm.
        % For the apical cell, the curve-fit parameters are as follows:
        %   imax = 2195.2 pA, x1 = 82.1 nm, α = 34.4 nm.

        x1 = (char_position(18e3) + char_position(13e3))/2;
        x2 = (char_position(1e3) + char_position(0.6e3))/2;

        % we assume linear dependence of the boltzmann parameters with
        % position perform linear regression
        % y = b0 + b1 * x
        x_jia2007 = [x1; x2];
        X = [ones(length(x_jia2007),1), x_jia2007];

        % boltzmann half-activatiobn parameter
        y_h = [35.3; 82.1];
        [h0, h1] = linear_regression(X, y_h);

        % boltzmann slope parameter
        y_s = [16.1; 34.4];
        [s0, s1] = linear_regression(X, y_s);

        MET_popen_BM_v1 = @(x,y) boltzmann(y, 1, h1*x + h0, s1*x + s0);
        MET_popen_BM_der_v1 = @(x,y) boltzmann_der(y, 1, h1*x + h0, s1*x + s0);
        
        if args.plotflag == true
            if isempty(args.handles) || numel(fieldnames(args.handles)) == 0
                hfig = figure();
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on
            
            u = linspace(-100e-9, 200e-9, 200);

            plot(u*1e9, MET_popen_BM_v1(zeros(size(u)),u*1e9), 'DisplayName', args.version);
            plot(u*1e9, MET_popen_BM_v1(ones(size(u)),u*1e9), 'DisplayName', args.version);
            xlabel('displacement (nm)')
            ylabel('open probability')
        end

        % -----------------------------------------------------------------
        
        % Johnson 2011 pp 1149
        % As with OHCs, perfusing 0.02 mM Ca2+ increased the peak size of
        % the MT current and also the fraction activated at rest (Figures 
        % 8A and 8B). The mean MT current increased from
        %    0.79 ± 0.07 nA (1.3 Ca2+; n = 4)
        % to 1.72 ± 0.12 nA (0.02 Ca2+; n = 5, T = 23?C)
        % and the fraction on at rest increased from 
        %    0.045± 0.004 (1.3 Ca2+) 
        % to 0.17 ± 0.03 (0.02 Ca2+).
        % Using the latter fraction and correcting the standing current to
        % 36?C yields a resting MTconductance of 5.7 nS.
        
        % My reasoning:
        % We want to shift boltzmann to be at rest (y=0) at p=0.17.
        % We have:
        %     p = (1 + exp(-(y-H)/S))^-1
        % then
        %     H = y + S ln(1/p - 1)
        
        p0 = 0.17;
        
        % H = @(x) 0 + (s1*x + s0) * log(1/p0 - 1);        
        
        % which transforms to        
        L = log(1./p0 - 1);
        H = @(x) (s1*x + s0) * L;
        
        MET_popen_BM_v2 = @(x,y) boltzmann(y, 1, H(x), s1*x + s0);
        MET_popen_BM_der_v2 = @(x,y) boltzmann_der(y, 1, H(x), s1*x + s0);
        
        MET_popen_BM_mag_v2 = @(x,y,m) boltzmann(y, m, H(x), s1*x + s0);
        MET_popen_BM_mag_der_v2 = @(x,y,m) boltzmann_der(y, m, H(x), s1*x + s0);
        
        HH = H(xgrid);
        SS = s1*xgrid + s0;
        
        MET_popen_BM_mag_pre_v2 = @(y,m) boltzmann(y, m, HH, SS);
        MET_popen_BM_mag_der_pre_v2 = @(y,m) boltzmann_der(y, m, HH, SS);
        
        if args.plotflag == true
            if isempty(args.handles) || numel(fieldnames(args.handles)) == 0
                hfig = figure();
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on
            
            u = linspace(-100e-9, 200e-9, 200);

            plot(u*1e9, MET_popen_BM_v2(zeros(size(u)),u*1e9), 'DisplayName', args.version);
            plot(u*1e9, MET_popen_BM_v2(ones(size(u)),u*1e9), 'DisplayName', args.version);
            xlabel('displacement (nm)')
            ylabel('open probability')
        end

        % -----------------------------------------------------------------
        
        % My reasoning:
        % We want to shift boltzmann to be at rest (y=0) at p=p(x).
        % We have:
        %     p = (1 + exp(-(y-H)/S))^-1
        % then
        %     H = y + S ln(1/p - 1)
        
        XX = [ones(2,1), [0; 1]];
        
        [q0, q1] = linear_regression(XX, [0.12; 0.22]);
        
        p0 = @(x) q1*x + q0;
        
        % H = @(x) 0 + (s1*x + s0) * log(1/p0 - 1);        
        
        % which transforms to        
        L = @(x) log(1./p0(x) - 1);
        H = @(x) (s1*x + s0) .* L(x);
        
        MET_popen_BM_v3 = @(x,y) boltzmann(y, 1, H(x), s1*x + s0);
        MET_popen_BM_der_v3 = @(x,y) boltzmann_der(y, 1, H(x), s1*x + s0);
        
        MET_popen_BM_mag_v3 = @(x,y,m) boltzmann(y, m, H(x), s1*x + s0);
        MET_popen_BM_mag_der_v3 = @(x,y,m) boltzmann_der(y, m, H(x), s1*x + s0);
        
        HH = H(xgrid);
        SS = s1*xgrid + s0;
        
        MET_popen_BM_mag_pre_v3 = @(y,m) boltzmann(y, m, HH, SS);
        MET_popen_BM_mag_der_pre_v3 = @(y,m) boltzmann_der(y, m, HH, SS);
        
        Numstacks = numel(xgrid);
        num_met_channels = 50;
        
        MET_popen_BM_mag_pre_v3_rand = @(y,m) m .* sum(rand(Numstacks,num_met_channels) < boltzmann(y, 1, HH, SS), 2) / num_met_channels;
        
        if args.plotflag == true
            if isempty(args.handles) || numel(fieldnames(args.handles)) == 0
                hfig = figure();
            else
                hfig = args.handles.popen;
                figure(hfig)
            end
            hold on
            
            u = linspace(-100e-9, 200e-9, 200);

            plot(u*1e9, MET_popen_BM_v3(zeros(size(u)),u*1e9), 'DisplayName', args.version);
            plot(u*1e9, MET_popen_BM_v3(ones(size(u)),u*1e9), 'DisplayName', args.version);
            xlabel('displacement (nm)')
            ylabel('open probability')
        end

        % -----------------------------------------------------------------
        
        ver = 'v3';
        ver = 'v2';
        
        switch ver
            case 'v1'
                met = @(varargin) MET_popen_BM_v1(varargin{:});
                met_der = @(varargin) MET_popen_BM_der_v1(varargin{:});
            case 'v2'
                met = @(varargin) MET_popen_BM_v2(varargin{:});
                met_der = @(varargin) MET_popen_BM_der_v2(varargin{:});
            case 'v3'
                met = @(varargin) MET_popen_BM_v3(varargin{:});
                met_der = @(varargin) MET_popen_BM_der_v3(varargin{:});
        end
        
        if args.plotflag
            hfig = figure;
            hold on
            xx = [0,1];

            plot(x_jia2007, y_h, 'ro', 'DisplayName', 'half activation (Jia 2007)');
            plot(xx, h1*xx + h0, 'r-', 'DisplayName', 'half activation  fit');

            plot(x_jia2007, y_s, 'bo', 'DisplayName', 'slope (Jia 2007)');
            plot(xx, s1*xx + s0, 'b-', 'DisplayName', 'slope fit');

            title('IHC boltzmann parameters')
            legend('location', 'NorthWest')
            
            if args.save_figures
                mySaveFig(hfig, 'MET_IHC_parameters_fit', [], 'Results/HairCells/', default_tikz_options{:});
            end
        end
        
        if args.plotflag
            hfig = figure;
            hold on
            xx = linspace(0,1,6);
            yy = linspace(-100, 200, 256);
            
            cmp = flip(rwb(numel(xx)));
            for i = 1:numel(xx)
                
                plot(yy, met(xx(i), yy), ...
                    'Color', cmp(i,:), ...
                    'LineStyle', '-', ...
                    'DisplayName', sprintf('%s', Frequency(char_frequency(xx(i)), 'Hz')));
            end
            % legend('location', 'SouthEast');

            p1 = plot(yy, boltzmann(yy, 1, y_h(1), y_s(1)), 'Color', 0.2 * ones(1,3), ...
                'LineStyle', '-.', 'DisplayName', 'Jia (2007) 15.5 kHz');
            p2 = plot(yy, boltzmann(yy, 1, y_h(2), y_s(2)), 'Color', 0.2 * ones(1,3), ...
                'LineStyle', '--', 'DisplayName', 'Jia (2007) 800 Hz');
            
%             plot([yy(1), 0], [p0, p0], ...
%                 'Color', 'k', ...[0.7, 0.7 0.7], ...
%                 'LineStyle', ':', ...
%                 'LineWidth', 0.75);
            
            legend([p1,p2], 'location', 'SouthEast')
            
            xlabel('BM displacement [nm]')
            ylabel('p_{open}')
            title('IHC MET open probability')
            
            cbarTL = cellfun(@(x) sprintf('%d', round(Frequency(char_frequency(x), 'Hz').Hz)), num2cell(xx), 'UniformOutput', false);
            
            manual_colorbar(cmp, cbarTL, ...
                'orientation', 'vertical', ...
                'reversed', true, ...
                'label', sprintf('f (Hz)'));
            
            if args.save_figures
                mySaveFig(hfig, 'MET_IHC', [], 'Results/HairCells/', default_tikz_options{:});
            end
            
        end
        
        if args.plotflag
            hfig = figure;
            hold on
            xx = linspace(0,1,6);
            yy = linspace(-100, 200, 256);
            
            cmp = flip(rwb(numel(xx)));
            for i = 1:numel(xx)
                
                plot(yy, met_der(xx(i),yy), ...
                    'Color', cmp(i,:), ...
                    'LineStyle', '-', ...
                    'DisplayName', sprintf('%s', Frequency(char_frequency(xx(i)), 'Hz')));
            end
            % legend('location', 'SouthEast');
            
            xlabel('BM displacement [nm]')
            ylabel('dp_{open}/d\xi')
            title('IHC MET open probability derivative')
            
            cbarTL = cellfun(@(x) sprintf('%d', round(Frequency(char_frequency(x), 'Hz').Hz)), num2cell(xx), 'UniformOutput', false);
            
            manual_colorbar(cmp, cbarTL, ...
                'orientation', 'vertical', ...
                'reversed', true, ...
                'label', sprintf('f (Hz)'));
            
            if args.save_figures
                mySaveFig(hfig, 'MET_IHC_derivative', [], 'Results/HairCells/', default_tikz_options{:});
            end
        end
        
        % -----------------------------------------------------------------
        
        switch ver
            case 'v2'
                MET_popen_BM_mag = MET_popen_BM_mag_v2;
                MET_popen_BM_mag_pre = MET_popen_BM_mag_pre_v2;

                MET_popen_BM_mag_der = MET_popen_BM_mag_der_v2;
                MET_popen_BM_mag_der_pre = MET_popen_BM_mag_der_pre_v2;

                MET_popen_BM = MET_popen_BM_v2;
                MET_popen_BM_der = MET_popen_BM_der_v2;
            case 'v3'
                MET_popen_BM_mag = MET_popen_BM_mag_v3;
                MET_popen_BM_mag_pre = MET_popen_BM_mag_pre_v3;
                % MET_popen_BM_mag_pre = MET_popen_BM_mag_pre_v3_rand;

                MET_popen_BM_mag_der = MET_popen_BM_mag_der_v3;
                MET_popen_BM_mag_der_pre = MET_popen_BM_mag_der_pre_v3;

                MET_popen_BM = MET_popen_BM_v3;
                MET_popen_BM_der = MET_popen_BM_der_v3;
        end
        
        % -----------------------------------------------------------------
        % IHC MET conductance
        
        % Maximal MET conductance is 29.8 ± 0.26 and 28.9 ± 0.65 nS 
        % for the apical- and basal-turn IHCs, respectively.
        
        % maximal conductance parameter
        y = [28.9; 29.8];

        % y = y * 1.6765; % to achieve resting conductance of 5.7 nS
        % y = y * 1.1400; % to achieve resting conductance of 5.7 nS

        [g0, g1] = linear_regression(X, y);
        
        MET_conductance_max = @(x) g1*x + g0;

        G_max_pre = MET_conductance_max(xgrid);
        
        FAC = 1;

        MET_conductance = @(x,y) MET_popen_BM_mag(x, FAC * y, g1*x + g0);
        MET_conductance_pre = @(y) MET_popen_BM_mag_pre(FAC * y, G_max_pre);
        MET_conductance_der = @(x,y) FAC * MET_popen_BM_mag_der(x, FAC * y, g1*x + g0);        
        MET_conductance_der_pre = @(y) FAC * MET_popen_BM_mag_der_pre(FAC * y, G_max_pre);
        
        if args.plotflag
            hfig = figure;
            hold on
            plot(xgrid, MET_conductance_max(xgrid), 'DisplayName', 'maximal conductance')
            plot(xgrid, MET_conductance(xgrid, zeros(size(xgrid))), 'DisplayName', 'resting conductance')
            title('IHC MET conductance')
            legend()
            xlabel('relative distance from stapes')
            ylabel('conductance (nS)')
        end
        
        
    case 'Johnson_2011'
        % CHECK THSESE VALUES (probably tweaked)

        % Johnson 2011 (pp. 1149)
        MET_conductance_rest = 5.7; % nS %

        % G_met_rest = @(x) linf(1.9, 3.6, x);
        % GAI_max = @(x) G_met_rest(x) / 0.18;

        MET_conductance_max = @(x) linf(10.5, 18, x);
        
        
    otherwise
        error('Unknown version %s', args.version)
end

    function [b, a] = linear_regression(X, y)
        
        beta = X\y;
        
        a = beta(2);
        b = beta(1);
    end

end

