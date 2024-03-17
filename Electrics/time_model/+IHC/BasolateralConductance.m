function [ G, dG, G_pre, dG_pre, dep, dep_pre, channels, EK_gradient ] = BasolateralConductance( xgrid, simulation_units, args )
%BASOLATERALCONDUCTANCE
%IHC

arguments
    xgrid (:,1) double = linspace(0, 1, 300)'
    simulation_units (1,1) struct = struct('voltage', 'V')
    args.version = 'Johnson_2015'
    args.plotflag = false
    args.save_figures = false
end

%%

% [char_frequency, char_position] = frequency_map();
[char_position, char_frequency] = characteristic_frequency_model('Li');

default_tikz_options = { ...
    'relativeDataPath', 'img/HairCells/', ...
    'extraaxisoptions', {'legend style={font=\tiny}'}, ...
    'parseStrings', true};

sim_to_mV = Unit.conversionConstant(simulation_units.voltage, 'mV');
mV_to_sim = Unit.conversionConstant('mV', simulation_units.voltage);

%%

% Data from Johnson 2011, Figure 8.E
% Parameters for the Boltzmann function:
% Voltage dependence of the K+ channel conductance

channels = [];

EK_gradient = [];

switch args.version

    case 'Tichacek_2022_v1'
        load('channels_OT_v1.mat', 'channels_apex', 'channels_base', 'EK_gradient');

        channels = channels_apex;

        G_A = zeros(size(channels));
        G_B = zeros(size(channels));

        for i = 1:numel(channels)
            par = fieldnames(channels(i).parameters);
            for j = 1:numel(par)
                if strcmp(par{j}, 'G')
                    G_A(i) = channels_apex(i).parameters.G;
                    G_B(i) = channels_base(i).parameters.G;
                else
                    assert(channels_base(i).parameters.(par{j}) == channels_apex(i).parameters.(par{j}))
                end
            end
        end

        for i = 1:numel(channels)
            channels(i).parameters.x_grad.G = [G_B(i), G_A(i)];
            channels(i).parameters.x_grad.EK = EK_gradient;
        end


        G_A = G_A * 1e9;
        G_B = G_B * 1e9;

        R_A = 1 ./ G_A;
        R_B = 1 ./ G_B;

        G_pre = zeros(numel(xgrid), numel(channels));

        for i = 1:numel(channels)
            G_pre(:,i) = linf(G_B(i), G_A(i), xgrid);
            % G_pre(:,i) = 1 ./ linf(1/R_B(i), 1/R_A(i), xgrid);
        end

        for i = 1:numel(channels)
            conductance{i} = @(P, x) P .* linf(G_B(i), G_A(i), x); % nS
            % conductance{i} = @(P, x) P ./ linf(R_B(i), R_A(i), x); % nS
            conductance_der{i} = @(x) linf(G_B(i), G_A(i), x); % nS/dP
            % conductance_der{i} = @(x) 1 ./ linf(R_B(i), R_A(i), x); % nS/dP
    
            conductance_pre{i} = @(P) G_pre(:,i) .* P; % nS
    
            % TODO -- these functions are in fact constatnt - include in jac_0
            conductance_der_pre{i} = @(P) G_pre(:,i); % nS/dP
        end

        G = conductance;
        dG = conductance_der;

        G_pre = conductance_pre;
        dG_pre = conductance_der_pre;
        
        dep = repmat({{'popen', 'xpos'}}, 1, numel(channels));
        dep_pre = repmat({{'popen'}}, 1, numel(channels));
         
        EK_gradient = EK_gradient * 1e3;

    case 'LopezPoveda_2006_Johnson_2015'
        
        % Johnson 2015 measured basal (~30 kHz) and apical (~300 Hz) turn
        %  in gerbil IHCs.
        % hf ... high-freqency ... apical part of the cochlea
        % lf ... low-freqency ... basal part of the cochlea

        x_hf = char_position(30e3);
        % 30 kHz is outside range in human but it should give reasonable
        % fit anyway

        x_lf = char_position(300);
        
        % The overall slope conductance, measured at around the likely in 
        % vivo membrane potential (Ac: –54 mV; Bc: –64 mV: Figure 2D) from 
        % the steady- state current-voltage (I-V) curves, was found to be 
        % significantly larger in apical (56.3 ± 3.0 nS, n = 26, P<0.0001) 
        % than in basal IHCs (10.0 ± 0.4 nS, n = 27).
        
        M_hf = 300; % nS
        M_lf = 900; % nS
        
        ratio = M_lf / M_hf;
        
        x_data = [x_hf; x_lf];
        X = [ones(length(x_data),1), x_data];
        
        
        % Now use the channel dynamics from LopezPoveda 2006
        
        channels = BasolateralChannelsIHC(args.version);
        
        G1 = channels(1).parameters.G * 1e9; % in nS;
        
        % G_lf = G1;
        G_lf = G1 * ratio;
        G_hf = G1;
        % G_hf = G1 / ratio;
        
        y_G = [G_hf; G_lf];
        [G10, G11] = linear_regression(X, y_G);
        
        fG1 = @(x) G11*x + G10;
        
        G1_pre = fG1(xgrid);
        
        
        G2 = channels(2).parameters.G * 1e9; % in nS;

        % G_lf = G2;
        G_lf = G2 * ratio;
        G_hf = G2;
        % G_hf = G2 / ratio;
        
        y_G = [G_hf; G_lf];
        [G20, G21] = linear_regression(X, y_G);
        
        fG2 = @(x) G21*x + G20;
        
        G2_pre = fG2(xgrid);
        
        conductance = {@(P,x) fG1(x) .* P, @(P,x) fG2(x) .* P}; % nS
        conductance_der = {@(x) fG1(x), @(x) fG2(x)}; % nS/dP

        conductance_pre = {@(P) G1_pre .* P, @(P) G2_pre .* P}; % nS

        % TODO -- these functions are in fact constatnt - include in jac_0
        conductance_der_pre = {@(P) G1_pre, @(P) G2_pre}; % nS/dP 
        
        G = conductance;
        dG = conductance_der;

        G_pre = conductance_pre;
        dG_pre = conductance_der_pre;
    
        dep = repmat({{'popen', 'xpos'}}, 1, numel(channels));
        dep_pre = repmat({{'popen'}}, 1, numel(channels));
        
        if args.plotflag == true
            
            hfig = figure;
            hold on
            V = linspace(-120e-3, 0e-3, numel(xgrid));
            
            sum_G = zeros(numel(V), 1);
            
            for j = 1:numel(channels)
                [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(j));
                
                single_G = Oinf .* channels(j).parameters.G * 1e9;
                sum_G = sum_G + single_G;
                
                plot(V*1e3, Oinf .* channels(j).parameters.G * 1e9, ...
                    'DisplayName', channels(j).name);
            end
            
            plot(V*1e3, sum_G, 'DisplayName', 'overall conductance');
            
            title('IHC basolateral conductance');
            xlabel('Voltage (mV)')
            ylabel('Conductance (nS)')
            legend('Location', 'NorthWest');
            
        end

        
    
    case {'LopezPoveda_2006', 'Dierich_2020'}
        
        channels = BasolateralChannelsIHC(args.version);
        
        switch args.version
            case 'LopezPoveda_2006'
                G1 = channels(1).parameters.G * 1e9; % in nS;
                G2 = channels(2).parameters.G * 1e9; % in nS;

                conductance = {@(P) G1 * P, @(P) G2 * P}; % nS
                conductance_der = {G1, G2}; % nS/dP
                
            case 'Dierich_2020'
                G1 = channels(1).parameters.G * 1e9; % in nS;
                G2 = channels(2).parameters.G * 1e9; % in nS;
                G3 = channels(3).parameters.G * 1e9; % in nS;
                G4 = channels(4).parameters.G * 1e9; % in nS;
                G5 = channels(5).parameters.G * 1e9; % in nS;
                G6 = channels(6).parameters.G * 1e9; % in nS;

                conductance = {@(P) G1 * P, @(P) G2 * P, @(P) G3 * P, @(P) G4 * P, @(P) G5 * P, @(P) G6 * P}; % nS
                conductance_der = {G1, G2, G3, G4, G5, G6}; % nS/dP
        end
        
        G = conductance;
        dG = conductance_der;

        G_pre = conductance;
        dG_pre = conductance_der;
    
        dep = repmat({{'popen'}}, 1, numel(channels));
        dep_pre = repmat({{'popen'}}, 1, numel(channels));
        
        if args.plotflag == true
            
            hfig = figure;
            hold on
            V = linspace(-120e-3, 0e-3, 200);
            
            sum_G = zeros(numel(V), 1);
            
            for j = 1:numel(channels)
                [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(j));
                
                single_G = Oinf .* channels(j).parameters.G * 1e9;
                sum_G = sum_G + single_G;
                
                plot(V*1e3, Oinf .* channels(j).parameters.G * 1e9, ...
                    'DisplayName', channels(j).name);
            end
            
            plot(V*1e3, sum_G, 'DisplayName', 'overall conductance');
            
            title('IHC basolateral conductance');
            xlabel('Voltage (mV)')
            ylabel('Conductance (nS)')
            legend('Location', 'NorthWest');
            
        end
        
    
    case {'Johnson_2011_v1', 'Johnson_2011_v2'}
                       
        switch args.version
            case 'Johnson_2011_v1'
                M = 580; % nS,  saturation value
                H = -26; % mV,  position of the inflexion point
                S = 11; % mV,   slope
            case 'Johnson_2011_v2'
                % alternative version (see Johnson)
                M = 470; % nS
                H = -31; % mV
                S = 10.5; % mV
        end
        
        H = H * mV_to_sim;
        S = S * mV_to_sim;
        
        conductance = @(V) boltzmann(V, M, H, S);
        conductance_der = @(V) boltzmann_der(V, M, H, S);
        
        G = conductance;
        dG = conductance_der;

        G_pre = conductance;
        dG_pre = conductance_der;
    
        dep = {'voltage'};
        dep_pre = {'voltage'};
        
    case 'Johnson_2015'
        
        % Johnson 2015 measured basal (~30 kHz) and apical (~300 Hz) turn
        %  in gerbil IHCs.
        % hf ... high-freqency ... apical part of the cochlea
        % lf ... low-freqency ... basal part of the cochlea

        x_hf = 0; % 30 kHz is outside range in human
        x_lf = char_position(300);
        
        % The overall slope conductance, measured at around the likely in 
        % vivo membrane potential (Ac: –54 mV; Bc: –64 mV: Figure 2D) from 
        % the steady- state current-voltage (I-V) curves, was found to be 
        % significantly larger in apical (56.3 ± 3.0 nS, n = 26, P<0.0001) 
        % than in basal IHCs (10.0 ± 0.4 nS, n = 27).
        
        M_hf = 300; % nS
        M_lf = 900; % nS
        
        H = -26; % mV
        S = 10.5; % mV
        
        H = H * mV_to_sim;
        S = S * mV_to_sim;
        
        if args.plotflag
            figure
            hold on
            plot(-54.0, 56.3, 'rx')
            plot(-64.0, 10.0, 'bx')

            V = linspace(-100, 0, 256);

            % plot(V, boltzmann(V, 470, -31, 10.5), 'k-');
            % plot(V, boltzmann(V, 580, -26, 11), 'k--');
        
            plot(V, boltzmann(V * mV_to_sim, M_hf, H, S), 'b-');
            plot(V, boltzmann(V * mV_to_sim, M_lf, H, S), 'r-');
        end
        
        x_data = [x_hf; x_lf];
        X = [ones(length(x_data),1), x_data];
        
        y_M = [M_hf; M_lf];
        [M0, M1] = linear_regression(X, y_M);
        
        M = @(x) M1*x + M0;
                
        G = @(V, x) boltzmann(V, M(x), H, S);
        dG = @(V, x) boltzmann_der(V, M(x), H, S);
        
        M_pre = M(xgrid);

        G_pre = @(V) boltzmann(V, M_pre, H, S);
        dG_pre = @(V) boltzmann_der(V, M_pre, H, S);
        
        conductance = G;
        conductance_pre = G_pre;
        
        dep = {'voltage', 'xpos'};
        dep_pre = {'voltage'};
        
        if args.plotflag
            hfig = figure;
            hold on
            
            xx = linspace(0,1,6);
            V = linspace(-120e-3, 0e-3, 200);
            
            for i = 1:numel(xx)
                plot(V*1e3, conductance(V, xx(i)))
            end
            
            title('IHC basolateral conductance');
            ylabel('Conductance (nS)')
            xlabel('Voltage (mV)')

        end
        
    otherwise
        error('Unknown version %s', args.version)
end



% GBI = @(x,V) boltzmann(V * sim_to_mV, M, H, S);
% dGBI = @(x,V) boltzmann_der(V * sim_to_mV, M, H, S) / mV_to_sim;



%%

if false %args.plotflag
    limits = [-90, 20];
    V = linspace(limits(1), limits(2), 100);
    
    Vss = -55; % estimate of the steady-state IHC voltage
    
    hfig = figure;
    hold on
    plot(V, conductance(V * mV_to_sim))
    plot(Vss, conductance(Vss * mV_to_sim), 'o', 'Color', 'black');
    frac = 100*conductance(Vss * mV_to_sim)/conductance(inf);
    text(Vss, conductance(Vss * mV_to_sim), ...
        sprintf(' \\leftarrow P(%g mV) = %0.2g %%', Vss, frac), ...
        'HorizontalAlignment', 'left')

    title('IHC G_K');
    xlabel('V [mv]');
    ylabel('G [nS]');
    
    if args.save_figures
        mySaveFig(hfig, 'IHC_G_K_fit', [], 'Results/HairCells/', default_tikz_options{:});
    end
end

%%

if args.plotflag
    limits = [-90,20];
    V = linspace(limits(1), limits(2), 200);
    
    [IHC_V_rest, ~, ~] = IHC.RestingPotential(xgrid, simulation_units);

    
    hfig = figure;
    hold on
   
    xx = linspace(0,1,10);
    
    for k = 1:numel(conductance)
                    
        cmp = flip(rwb(numel(xx)));
        for i = 1:length(xx)            
            Vss = IHC_V_rest(xx(i));
            
            dep_eval = {};
            dep_eval_SS = {};
            for l = 1:numel(dep{k})
                switch dep{k}{l} 
                    case 'popen'
                        [~,~,~,~,~,~,~,dep_eval{end+1}] = pOpenIHC(V * mV_to_sim, channels(k), 'mass', false, 'system', false, 'rhs', true, 'jacobian', false);
                        [~,~,~,~,~,~,~,dep_eval_SS{end+1}] = pOpenIHC(Vss * mV_to_sim, channels(k), 'mass', false, 'system', false, 'rhs', true, 'jacobian', false);
                    case 'voltage'
                        dep_eval{end+1} = V * mV_to_sim;
                        dep_eval_SS{end+1} = Vss * mV_to_sim;
                    case 'xpos'
                        dep_eval{end+1} = xx(i);
                        dep_eval_SS{end+1} = xx(i);
                end
            end
            
            plot(V, conductance{k}(dep_eval{:}), 'Color', cmp(i,:), 'LineStyle', '-.', 'Marker', 'none')
            plot(Vss, conductance{k}(dep_eval_SS{:}), 'Color', cmp(i,:), 'LineStyle', 'none', 'Marker', 'o')
        end
    end
    
    plot(-64.0, 10.0, 'kx') % base - hf
    plot(-54.0, 56.3, 'kx') % apex - lf
    
    
    title('IHC basolateral conductance');
    xlabel('V [mv]');
    ylabel('IHC G_K [nS]');
    
    if false
        cbarTL = cellfun(@(x) sprintf('%d', round(Frequency(char_frequency(x), 'Hz').Hz)), num2cell(xx), 'UniformOutput', false);

        manual_colorbar(cmp, cbarTL, ...
            'orientation', 'vertical', ...
            'reversed', true, ...
            'label', sprintf('f (Hz)'));
    end
    
    if args.save_figures
        mySaveFig(hfig, 'IHC_G_K_fit', [], 'Results/HairCells/', default_tikz_options{:});
    end
end

%%

    function [b, a] = linear_regression(X, y)
        
        beta = X\y;
        
        a = beta(2);
        b = beta(1);
    end

end

