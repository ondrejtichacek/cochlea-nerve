function [ nonlinParams ] = amplifier_params( K, args )
%AMPLIFIER_PARAMS
arguments
    K (:,1) double
    args.plotflag (1,1) logical = false
    args.hfig = []
    args.y_range (1,2) double = [NaN, NaN];
    args.K_range (1,2) double = [NaN, NaN];
    args.X_plot_scale (1,1) double = 1
    args.Y_plot_scale (1,1) double = 1
end
if all(isnan(args.y_range))
    M = max(K);
    args.y_range = [-M, M];
end
if all(isnan(args.K_range))
    M = max(K);
    m = min(K);
    args.K_range = [m, M];
end

% Original values used in the models by NOBILI and VETESNIK:

% NOBILI
% nonlinParams = struct( ...
%     'y1', 0.1139, ...
%     'y2', 0.3736, ...
%     'c1', 0.7293, ...
%     'c2', 1.4974, ...
%     'b', 0.30991, ...
%     'q', 1 )

% VETESNIK
% nonlinParams = struct( ...
%     'y1', 0.01139, ...
%     'y2', 0.03736, ...
%     'c1', 0.7293, ...
%     'c2', 1.4974, ...
%     'b', 0.30991, ...
%     'q', 0.1 )

% The above can be obtained as follows:
% nonlinParams = amplifier_params(1); % VETESNIK
% nonlinParams = amplifier_params(0.1); % NOBILI

% NEW PARAMETERS
% nonlinParams = struct( ...
%     'y1', 0.113538137421216, ...
%     'y2', 0.375319620486122, ...
%     'c1', 0.728962829291601, ...
%     'c2', 1.497538031887178, ...
%     'b', 0.309933278741342, ...
%     'q', 0.8 )

% -------------------------------------------------------------------------

for i = 0:100
    if i == 0 || any(isnan(b))
        nonlinParams = struct( ...
            'y1', (K + i*eps) * 0.113538137421216, ...
            'y2', (K + i*eps) * 0.375319620486122, ...
            'c1', 0.728962829291601, ...
            'c2', 1.497538031887178, ...
            'b', 0.309933278741342, ...
            'q', (K + i*eps) * 1 ...
            );
    else
        break
    end
    
    % evaluate without b
    b = boltzmann_3level_wrapper(0, nonlinParams);
    
end
    
assert(~any(isnan(b)));

nonlinParams.b = b;

% now it should be zero

    function Y = boltzmann_3level_wrapper(y, nonlinParams)
        
        Y = boltzmann_3level(y, nonlinParams.c1, nonlinParams.c2, nonlinParams.y1, nonlinParams.y2);
    end

if args.plotflag == true
    if isempty(args.hfig)
        args.hfig = figure;
    else
        actfigure(args.hfig);
    end
    hold on
        
    k = linspace(args.K_range(1), args.K_range(2), 10);
    c = jet(10);
    
    y = linspace(args.y_range(1), args.y_range(2), 200);
    
    for i = 1:numel(k)
        nlp = amplifier_params(k(i));
    
        Y = nonlin(y, nlp.y1, nlp.y2, nlp.c1, nlp.c2, nlp.b, nlp.q);
        plot(args.X_plot_scale * y, args.Y_plot_scale * Y, 'Color', c(i,:));
    end
end

end
