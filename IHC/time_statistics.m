function [ stats ] = time_statistics( simulation_result, variables_to_analyse, runopt )
%TIME_STATISTICS 


time = simulation_result.t.ms;

% find time index of steady state
initTimeIndex = 1;

% find time index of steady state
steadyStateTime = runopt.zeroTime.ms;
steadyStateTimeIndex = min([length(time), find(steadyStateTime <= time,1)]);


% find first time index of analysis
analysisStartTime = runopt.analysisStart.ms;
analysisStartTimeIndex = min([length(time), find(analysisStartTime <= time,1)]);

% restrict variables to analysis time
time = time(analysisStartTimeIndex:end);


stats = variables_to_analyse;

variables_to_analyse = fields(stats);

for k = 1:numel(variables_to_analyse)
    variable_name = variables_to_analyse{k};
    
    % extract variable
    variable = simulation_result.(variable_name);
    
    % handle units
    if isa(variable, 'Unit')
        unit = variable.internal_unit;
        variable = variable.(unit);
    else
        unit = '';
    end    
    stats.(variable_name).unit = unit;
    
    stats.(variable_name).statistics = struct( ...
        'init', struct(), ...
        'steady_state', struct(), ...
        'tspan', struct());
    
    stats.(variable_name).statistics.init = ...
        analyse_signal(variable(initTimeIndex,:));

    stats.(variable_name).statistics.steady_state = ...
        analyse_signal(variable(steadyStateTimeIndex,:));
    
    stats.(variable_name).statistics.tspan = ...
        analyse_signal(variable(analysisStartTimeIndex:end,:));
    
end

    function res = analyse_signal(sig)
        res = struct;

        % first slice
        % S.subs = repmat({':'},1,ndims(sig));
        % S.subs{1} = 1; % the first row
        % S.type = '()';
        % res.val0 = subsasgn(sig,S,[]);
        res.val0 = sig(1);

        res.min = min(sig);
        res.max = max(sig);
        res.mean = mean(sig);
        res.std = std(sig);
    end

end

