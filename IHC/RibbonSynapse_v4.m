classdef RibbonSynapse_v4
    %RIBBONSYNAPSE_V4 
    
    properties
        rho
        psi
        
        xv
        xc
        method
        
        num_release_sites;
        num_channels;
        
        distance_vesicle_membrane
        intervesicle_distance
        
        channel_radius
        vesicle_radius
        
        channel_distribution_method
        channel_distribution_parameters

        area
    end
    
    methods (Static)
        % function k = TransmitterRelease(C, k_max, H, S)
        %     %TRANSMITTERRELEASE 
        % 
        %     k = boltzmann(C, k_max, H, S);
        % end
        function k = TransmitterRelease(C, k_max, n, KA)
            %TRANSMITTERRELEASE 
        
            k = k_max .* Hill_Langmuir_A(C, n, KA);
        end
    end
    
    methods
        function obj = RibbonSynapse_v4(args)
            %RIBBONSYNAPSE_V4
            % 74 Ca2+ channels (see text) and 14 docked vesicles (Wong et al., 2014).
            % The diameter of Ca2+ channels is ∼15 nm (Wolf et al., 2003), while that of vesicles is
            % ∼40 nm (Khimich et al., 2005). Placing 74 Ca2+ channels closely packed
            % (7.5 nm minimum distance) gives a presynaptic area of ∼0.033 µm2 , which is
            % consistent with the most recent EM estimates (0.066 µm2 ; Khimich et al.,
            % 2005) and STED microscopy (0.11 µm2 and 0.097 µm2 for apical and middle
            % cochlea IHCs, respectively: Meyer e al., 2009; 0.034 µm2 for mature IHCs:
            % Wong et al., 2014)
            arguments
                
                args.num_release_sites (1,1) double {mustBeInteger, mustBePositive} = 14
                args.num_channels (1,1) double {mustBeInteger, mustBePositive} = 5
                
                args.distance_vesicle_membrane (1,1) double {mustBePositive} = 20 % nm
                args.intervesicle_distance (1,2) double {mustBePositive} = [55,55] % nm

                args.channel_radius (1,1) double {mustBePositive} = 7.5 % nm
                args.vesicle_radius (1,1) double {mustBePositive} = 20 % nm
                
                args.channel_distribution_method (1,:) char = 'uniform_ring'
                
                args.channel_distribution_parameters (1,:) struct = ...
                    struct('dist_min', 7.5 + 15, ...
                           'dist_max', 40, ...
                           'exclusion_dist_vesicle', 7.5 + 15, ...
                           'exclusion_dist_channel', 2*7.5 + 1);
                       
                args.plotflag (1,1) logical = false
                args.plot_histograms (1,1) logical = true
            end
            
            f = fieldnames(args);
            f = setdiff(f, {'plotflag', 'plot_histograms'});
            for i = 1:numel(f)
                obj.(f{i}) = args.(f{i});
            end
            
            %%
            
            n = obj.num_release_sites;
            M = obj.num_channels;
            
            assert(mod(n, 2) == 0)
            
            obj.xv = NaN(n, 2);
            
            for i = 1:n
                obj.xv(i, 1) = floor((i-1)/2);
                obj.xv(i, 2) = mod(i,2);
            end
            
            obj.xv = obj.xv .* obj.intervesicle_distance;
            
            %% Channel distribution
            

            obj.xc = obj.channel_distribution(n, M);

            
            %% Add 3rd dimension
            
            obj.xv(:,3) = obj.distance_vesicle_membrane + obj.vesicle_radius;
            obj.xc(:,3) = 0;
            
            %% Compute distances
            
            % dist channel -- vesicle membrane
            obj.rho = zeros(n, M);
            for i = 1:n
                obj.rho(i,:) = sqrt(sum((obj.xc - obj.xv(i,:)).^2, 2)) - obj.vesicle_radius;
            end

            assert(all(obj.rho(:) >= 0));
            
            % dist channel -- channel
            obj.psi = zeros(M, M);
            for i = 1:M
                obj.psi(i,:) = sqrt(sum((obj.xc - obj.xc(i,:)).^2, 2));
            end
            
            a1 = min(obj.xc(:,1));
            b1 = max(obj.xc(:,1));
            
            a2 = min(obj.xc(:,2));
            b2 = max(obj.xc(:,2));
            
            area = (b1-a1)*(b2-a2);
            
            obj.area = area;
            
            %%
            
            if args.plotflag == true
            
                % figure
                % hold on
                % plot(obj.xv(:,1), obj.xv(:,2), 'bo')
                % plot(obj.xc(:,1), obj.xc(:,2), 'rx')
                % axis equal
                % 
                % figure
                % hold on
                % plot3(obj.xv(:,1), obj.xv(:,2), (0.4 + 0.2)*ones(n,1), 'bo')
                % plot3(obj.xc(:,1), obj.xc(:,2), zeros(M,1), 'rx')
                % axis equal

                figure
                hold on
                [X,Y,Z] = sphere;

                X = obj.vesicle_radius * X;
                Y = obj.vesicle_radius * Y;
                Z = obj.vesicle_radius * Z;

                for i = 1:n

                    if rand(1) < 4/14
                        c = 'b';
                    else
                        c = [0.3, 0.3, 0.3];
                    end

                    surf(X + obj.xv(i,1), Y + obj.xv(i,2), Z + obj.xv(i,3), ...
                        'FaceColor', c, ...
                        'EdgeColor', c, ...
                        'FaceAlpha',0.05, ...
                        'EdgeAlpha',0.15);
                end
                [xp, yp] = circle(obj.vesicle_radius);

                for i = 1:n
                    plot(xp + obj.xv(i,1), yp + obj.xv(i,2), ...
                        'Color', [.7,.7,.7])
                end

                % viscircles(obj.xv(:,1:2), obj.vesicle_radius*ones(1,n), ...
                %     'Color', [0,0,0.2])

                % plot3(obj.xc(:,1), obj.xc(:,2), zeros(M,1), 'rx')
                % viscircles(obj.xc(:,1:2), obj.channel_radius*ones(1,M), 'Color', 'r')

                [X,Y,Z] = cylinder;

                X = obj.channel_radius * X;
                Y = obj.channel_radius * Y;
                Z = 5 * Z;

                for i = 1:M
                    surf(X + obj.xc(i,1), Y + obj.xc(i,2), Z + obj.xc(i,3), ...
                        'FaceColor', 'r', ...
                        'EdgeColor', 'r', ...
                        'FaceAlpha',0.2, ...
                        'EdgeAlpha',0.4);
                end

                view(-24,36)
                zlim([0,100])
                grid on
                axis equal

                title('Ribbon Synapse')
                subtitle(sprintf('presynaptic membrane area = %.2g um^2\n area per channel = %.2g um^2\nconductance = %g nS (%.2g nS/um^2)', ...
                    area*1e-6, area/M*1e-6, 15*1e-3*M, 15*1e-3*M/(area*1e-6)))

                if args.plot_histograms
                    axes('Position',[.7 .05 .2 .2])
                    box on
                    h = histogram(obj.rho, 0:20:(max(obj.rho(:))+20));
                    h.FaceColor = [0.1 0 0.9];
                    xlabel('distance (nm)')
                    ylabel('count')
                    title('dist channel -- vesicle membrane')
    
                    axes('Position',[.45 .05 .2 .2])
                    box on
                    h = histogram(obj.psi(obj.psi>0), 0:20:(max(obj.psi(:))+20));
                    h.FaceColor = [0.9 0 0.1];
                    xlabel('distance (nm)')
                    ylabel('count')
                    title('dist channel -- channel')
                end
                drawnow
            
            end
        end
        function xc = channel_distribution(obj, n, M)
            
            param = obj.channel_distribution_parameters;

            xc = NaN(M, 2);
            
            switch obj.channel_distribution_method

                case 'LJ'
                    o = [mean([max(obj.xv(:,1)), min(obj.xv(:,1))]), ...
                         mean([max(obj.xv(:,2)), min(obj.xv(:,2))])];
                    
                    e = [[max(obj.xv(:,1)), min(obj.xv(:,1))]; ...
                         [max(obj.xv(:,2)), min(obj.xv(:,2))]];

                    s = param.s;

                    w = -diff(e(1,:)) + 2*s;
                    l = -diff(e(2,:)) + 2*s;

                    for k = 1:M    
                        x1 = l*(rand(1,1) - 0.5);
                        x2 = w*(rand(1,1) - 0.5);

                        X = o + [x1, x2];

                        xc(k,:) = X;
                    end

                    % vesicles
                    % epsilon0 = 1;
                    % r0 = 150;
                    % n0 = 4;
                    epsilon0 = param.epsilon0;
                    r0 = param.r0;
                    n0 = param.n0;
                    
                    % channels
                    % m1 = 10;
                    % h1 = 20;
                    % s1 = 4;
                    m1 = param.m1;
                    h1 = param.h1;
                    s1 = param.s1;

                    % epsilon1 = 0.05;
                    % r1 = 20;
                    % n1 = 6;

                    r_excl = 10; % nm

                    beta_jump = 10;

                    for it = 1:2000
                        for k = 1:M

                            X = xc(k,:);

                            if it < 10000
                                x1 = l*(rand(1,1) - 0.5);
                                x2 = w*(rand(1,1) - 0.5);    
                                Xnew = o + [x1, x2];
                            else
                                d = beta_jump * rand(1);
                                phi = 2*pi*rand(1);
                                jump = d * [sin(phi), cos(phi)];
                                % jump = beta_jump * (rand(1,2) - 0.5);
                                Xnew = X + jump;
                            end

                            xcc = xc;
                            xcc(k,:) = [];

                            phi_loc = sqrt(sum((obj.xv - X).^2, 2));
                            phi_loc_new = sqrt(sum((obj.xv - Xnew).^2, 2));

                            psi_loc = sqrt(sum((xcc - X).^2, 2));
                            psi_loc_new = sqrt(sum((xcc - Xnew).^2, 2));
                              
                            q = r_excl./psi_loc;
                            q(q > 1) = 1;

                            q_new = r_excl./psi_loc_new;
                            q_new(q_new > 1) = 1;

                            F = [ ...
                                lennard_jones(phi_loc, epsilon0, r0, n0);
                                phi_loc.^0.15 / n
                                q.^12
                                - boltzmann(psi_loc, m1, h1, s1);
                                ...boltzmann(phi_loc, 0.1, 100, 20);                                
                                ...14*lennard_jones(min(phi_loc), epsilon0, r0, n0);
                                ...- boltzmann(phi_loc, epsilon0, r0, n0);
                                ...lennard_jones(psi_loc, epsilon1, r1, n1);
                                ];
                            F_new = [ ...
                                lennard_jones(phi_loc_new, epsilon0, r0, n0);
                                phi_loc_new.^0.15 / n
                                q_new.^12
                                - boltzmann(psi_loc_new, m1, h1, s1);
                                ...14*lennard_jones(min(phi_loc_new), epsilon0, r0, n0);
                                ...lennard_jones(psi_loc_new, epsilon1, r1, n1);
                                ];

                            deltaE = sum(F_new) - sum(F);
                            % deltaE is negative if X_new better than X

                            beta = 50;
                            beta = 100;

                            if rand(1) < exp(-beta * deltaE)
                                xc(k,:) = Xnew;
                            end
                        end
                    end

                case 'uniform'
            
                    a = -0.3;
                    b = 0.3;
                    
                    X = x + a + (b-a)*rand(1,2);
                    
                case 'uniform_ring'
                    
                    a = param.dist_min;
                    b = param.dist_max;
                    
                    dc = param.exclusion_dist_channel;
                    dv = param.exclusion_dist_vesicle;
                    
                    m = M / n;
                    mustBeInteger(m);
                    
                    for i = 1:n
                        for j = 1:m
                            k = m*(i-1) + j;
                            
                            x = obj.xv(i,:);
                            
                            while true
                                dist = a + (b-a)*rand(1,1);
                                angle = 2*pi*rand(1,1);

                                X = x + dist * [cos(angle), sin(angle)];

                                psi_loc = sqrt(sum((xc - X).^2, 2));
                                phi_loc = sqrt(sum((obj.xv - X).^2, 2));

                                if ~any(psi_loc < dc) && ~any(phi_loc < dv)
                                    break
                                end

                            end
                            
                            xc(k,:) = X;
                        end
                    end
                    
                    
                case 'uniform_band'
                    
                    o = [mean([max(obj.xv(:,1)), min(obj.xv(:,1))]), ...
                         mean([max(obj.xv(:,2)), min(obj.xv(:,2))])];
                    
                    w = param.width;
                    l = param.length;
                    
                    dc = param.exclusion_dist_channel;
                    dv = param.exclusion_dist_vesicle;
                    
                    maxtries = 1000;

                    for k = 1:M    
                        i = 0;
                        while true
                            i = i + 1;
                            x1 = l*(rand(1,1) - 0.5);
                            x2 = w*(rand(1,1) - 0.5);

                            X = o + [x1, x2];

                            psi_loc = sqrt(sum((xc - X).^2, 2));
                            phi_loc = sqrt(sum((obj.xv - X).^2, 2));

                            if ~any(psi_loc < dc) && ~any(phi_loc < dv)
                                break
                            else
                                if i > maxtries
                                    error('Channel distribution not possible')
                                end
                            end

                        end

                        xc(k,:) = X;
                    end
                
                case 'regular_band'
                    
                    o = [mean([max(obj.xv(:,1)), min(obj.xv(:,1))]), ...
                         mean([max(obj.xv(:,2)), min(obj.xv(:,2))])];
                    
                    w = param.width;
                    l = param.length;
                    
                    x = linspace(-l/2, l/2, ceil(M/3));
                    y = linspace(-w/2, w/2, 3);

                    [X, Y] = meshgrid(x,y);

                    xc = [X(:), Y(:)];                        
                    xc = o + xc(1:M,:);
                
                case 'regular_band_4'
                    
                    N = 4;

                    o = [mean([max(obj.xv(:,1)), min(obj.xv(:,1))]), ...
                         mean([max(obj.xv(:,2)), min(obj.xv(:,2))])];
                    
                    w = param.width;
                    l = param.length;
                    
                    x = linspace(-l/2, l/2, ceil(M/N));
                    y = linspace(-w/2, w/2, N);

                    [X, Y] = meshgrid(x,y);

                    xc = [X(:), Y(:)];
                    xc = o + xc(1:M,:);
                
                case 'regular_band_n'
                    
                    N = numel(param.wpos);

                    o = [mean([max(obj.xv(:,1)), min(obj.xv(:,1))]), ...
                         mean([max(obj.xv(:,2)), min(obj.xv(:,2))])];
                    
                    l = param.length;
                    
                    x = linspace(-l/2, l/2, ceil(M/N));
                    y = param.wpos;

                    [X, Y] = meshgrid(x,y);

                    xc = [X(:), Y(:)];
                    xc = o + xc(1:M,:);

                otherwise
                    error('Not set up')
            end
        end
    end
end

