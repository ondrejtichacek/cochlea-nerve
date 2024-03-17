function H = plot_steady_state_and_profile( xx, Y, analysisStartFrame, SCALING, relative, varargin )
%PLOT_STEADY_STATE_AND_PROFILE 
arguments
    xx
    Y
    analysisStartFrame (1,1) double = 1
    SCALING (1,1) double = 1
    relative (1,1) logical = false
end
arguments (Repeating)
    varargin
end

Y0 = Y(1,:);
Ymean = mean(Y(analysisStartFrame:end,:),1);

if relative == true
    SS = Y0;
else
    SS = Ymean;
end

% ERR = std(Y - SS, 0, 1);

O = zeros(size(SS));

ERR = [ ...
    abs(max(max(Y(analysisStartFrame:end,:) - SS, [], 1), O)); ...
    abs(min(min(Y(analysisStartFrame:end,:) - SS, [], 1), O))];

if SCALING ~= 1
    SS = SS * SCALING;
    ERR = ERR * SCALING;
end

if relative == true
    SS = O;
end

H = shadedErrorBar(xx, SS, ERR, 'transparent', true, 'lineProps', 'k:', varargin{:});

hold on

% plot steady state
plot(xx, Y0, 'k-');

end

