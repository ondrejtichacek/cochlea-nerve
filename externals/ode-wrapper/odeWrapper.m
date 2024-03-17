function [ varargout ] = odeWrapper( wrapopt, SOLVER, ode, tspan, y0, options, varargin )
%ODEWRAPPER
%
% IS THIS MATHEMATICALLY CORRECT???

if ischar(SOLVER)
    SOLVER = str2func(SOLVER);
end

if isempty(wrapopt)
    % default (unchanged call)
    [varargout{1}, varargout{2}] = SOLVER(ode, tspan, y0, options, varargin{:});
    
elseif isfield(wrapopt, 'v2') && wrapopt.v2 == true
    
    if isfield(wrapopt, 'InterDivFcn')
        InterDivVars = wrapopt.InterDivFcn(tspan);
    else
        InterDivVars = {};
    end
    
    timer = tic;
        
    SOLVER(ode, tspan, y0, options, varargin{:}, InterDivVars{:});
    
    statistics = wrapopt.storage.statistics;
    
    statistics.global_time = toc(timer);
    statistics.ode_time = statistics.global_time - statistics.transfer_time - statistics.save_time; 
    
    varargout{1} = statistics;
    
else
    
    switch wrapopt.save_method
        case {'matlab_matfile', 'struct'}
            ODEWRITE = @matlab_matfile;
        case 'matlab_fwrite'
            ODEWRITE = @matlab_fwrite;
        case 'c_posix'
            ODEWRITE = @c_posix;
    end
    
    assert(~isempty(wrapopt.NumDiv))
    % Divide the integration to `NumDiv` steps to lower memory requirements
    
    
    if ischar(wrapopt.NumDiv) && strcmp(wrapopt.NumDiv, 'auto')
        if numel(tspan) <= 2
            error('wrapopt.NumDiv = auto option is not available for unspecified time samples');
        end
        
        memReq = numel(y0) * sizeof(class(y0)) * numel(tspan);
        memLim = wrapopt.memopt.maxVarMem;
        
        wrapopt.NumDiv = ceil(memReq/memLim);
        
    end

    % If time samples are not provided, but options.TimeStep is, create the
    % time samples manually
    if numel(tspan) == 2
        if isa(options, 'solverOpt') && ~isempty(options.TimeStep)
            % create equidistant output samples
            numSteps = ceil((tspan(2) - tspan(1)) / options.TimeStep);
            numSamples = numSteps + 1;

            tspan = linspace(tspan(1), tspan(2), numSamples);
            tspan = tspan(:);
        end
    end
    
    % First create time intervals for consecutive integration:
    % Distinguish between the two cases
    %   tspan = [t0, tf]
    % and
    %   tspan = [t0, t1, t2, ..., tf]
    %
    if numel(tspan) == 2
        % In the first case (assuming t0 = 0 and NumDiv = 9)
        %   T = { [      0 ,   tf/9],
        %         [   tf/9 , 2*tf/9],
        %              ...
        %         [ 8*tf/9 , 2*tf/9] }
        %
        tt = linspace(tspan(1), tspan(2), wrapopt.NumDiv + 1);
        T{numel(tt)-1} = [];
        for i = 1:numel(tt)-1
            T{i} = [tt(i), tt(i+1)];
        end
    else
        % In the second case (assuming tspan = [0, 0.1, ..., 9] and NumDiv = 9)
        %   T = { [ 0, 0.1, 0.2, ..., 1 ],
        %         [ 1, 1.1, 1.2, ..., 2 ],
        %              ...
        %         [ 8, 8.1, 8.2, ..., 9 ], }
        %
        k = ceil((numel(tspan) - 1) / wrapopt.NumDiv);
        tt = unique([1:k:numel(tspan), numel(tspan)]);
        T{numel(tt)-1} = [];
        for i = 1:numel(tt)-1
            T{i} = tspan(tt(i) : tt(i+1));
        end
    end

    % In the second case, it may be the case that the desired NumDiv is
    % too high (larger than or equal to the number of elements of tspan).
    % In that case, presenet a warning.
    %
    if numel(T) ~= wrapopt.NumDiv
        warning('Desired number of divisions %d replaced by %d\n', wrapopt.NumDiv, numel(T));
        warning('!!! You should set this in the saved options !!!');
    end

    % Handle checkpoints
    if wrapopt.save_checkpoints
        
        if numel(tspan) == 2
            warning('Generated time samples are not necessarily equidistant, the checkpoint frequency/period may not be as desired')
        end
        
        % Get dt of time samples
        if isa(options, 'solverOpt') && ~isempty(options.TimeStep)
            dt = options.TimeStep;
        else
            dt = tspan(2) - tspan(1);
        end
        
        % Requested period of checkpoints in sense of indices
        checkpoint_period_ind = max(1, round(wrapopt.checkpoint_period / dt));
        
        % Checkpoints time indices for each interval
        checkpoint_indices{numel(T)} = [];
        for i = 1:numel(T)
            checkpoint_indices{i} = 1 : checkpoint_period_ind : numel(T{i});
            
            % Allways remove first index and add last
            checkpoint_indices{i} = unique([checkpoint_indices{i}(2:end), numel(T{i})]);
        end
    end
    
    % Now perform the actual integration steps:
    % We need to track the last index of yet written solution `nt0`.
    
    if wrapopt.append == true
        if exist(wrapopt.res.t.Properties.Source, 'file')
            switch wrapopt.save_method
                case 'matlab_matfile'
                    nt0 = max(matfileWhosByName(wrapopt.res.t, 't', 'size')); % = memmory friendly length(wrapopt.res.t.t)
                    nt0  = nt0 - 1;
                case 'matlab_fwrite'
                    error('Not yet implemented')
                case 'c_posix'
                    error('Not yet implemented')
            end
        else
            nt0 = 0;
        end
    else
        % Initialize its value to 0 as nothing have been written yet.
        nt0 = 0;
    end
    
    assert(size(y0,2) == 1, 'Vector of initial conditions must be a column vector');

    % write initial condition
    nt0 = ODEWRITE( tspan(1), y0', wrapopt.res, nt0);
    
    % Save initial condition as first checkpoint
    if wrapopt.save_checkpoints
        
        % time index of last checkpoint
        nt0c = wrapopt.checkpoints_last_t_index;
        
        % if nothing has been written to checkpoints
        if nt0c == 0
            nt0c = ODEWRITE( tspan(1), y0', wrapopt.checkpoints, nt0c);
        end
    end
    
    statistics.ode_time = NaN(1,numel(T));
    statistics.save_time = NaN(1,numel(T));

    for i = 1:numel(T) % numel(T) is the actual number of integration steps

        timerOde = tic;

        % If applicable, compute variables that depend on current tspan
        % e.g. load inputs for current tspan
        if isfield(wrapopt, 'InterDivFcn')
            InterDivVars = wrapopt.InterDivFcn(T{i});
        else
            InterDivVars = {};
        end
        
        % call the solver setting `tspan = T{i}`
        switch func2str(SOLVER)
            case {'BDF2Solver', 'odeEuler'}
                [t, y] = SOLVER(ode, T{i}, y0, options, varargin{:}, InterDivVars{:});
            otherwise
                if i > 1
                    options = odeset(options, 'InitialStep', diff(sol.x(end-1:end)));
                end
                sol = SOLVER(ode, T{i}, y0, options, varargin{:}, InterDivVars{:});
                
                t = T{i};
                y = deval(sol, t)';
                
        end
        
        statistics.ode_time(i) = toc(timerOde);
        
        timerSave = tic;
        
        % write the solution: `2:end` because the last time point of
        % previous iteration = first time point of current iteration
        nt0 = ODEWRITE( t(2:end), y(2:end,:), wrapopt.res, nt0);
        
        if wrapopt.save_checkpoints
            
            % get the indices of valid checkpoints (becasue if the solver
            % fails and t, y are incomplete, we would get a bad subscript
            % error)
            sel = checkpoint_indices{i} <= numel(t);
            ind = checkpoint_indices{i}(sel);
            
            % write checkpoints to file
            if ~isempty(ind)
                nt0c = ODEWRITE( t(ind), y(ind,:), wrapopt.checkpoints, nt0c);
            end
        end

        statistics.save_time(i) = toc(timerSave);
        
        % update the initial conditions for the next step
        y0 = y(end,:)';
        
        if wrapopt.verbose >= 2
            fprintf('The integration step %d/%d completed in %s.\n', ...
                i, numel(T), disp_toc(toc(timerOde)));
        end
        
        % check if the solver returned the whole solution (matlab sometimes
        % only throws a warning and returns incomplete solution)
        if abs(t(end) - T{i}(end)) > eps
                        
            % if there is a fallback solver specified, we will try to
            % switch to it for some time (100 steps) and then switch back
            % to the original solver
            if isfield(wrapopt, 'fallback_solver') && ~isempty(wrapopt.fallback_solver)
                warning('Failure occured during odeWrapper (solver %s), trying fallback solver (%s)', func2str(SOLVER), wrapopt.fallback_solver)
                
                i1 = numel(t);
                i2 = min(numel(T{i}),i1+100);
                
                % tspan for the fallback solver
                tspan_fallback = T{i}(i1:i2);
                
                % initial condition for the fallback solver
                y0_fallback = y(end,:)';
                
                SOLVER_FALLBACK = wrapopt.fallback_solver;
                if ischar(SOLVER_FALLBACK)
                    SOLVER_FALLBACK = str2func(SOLVER_FALLBACK);
                end
                
                % If applicable, compute variables that depend on current tspan
                % e.g. load inputs for current tspan
                if isfield(wrapopt, 'InterDivFcn')
                    InterDivVars_fallback = wrapopt.InterDivFcn(tspan_fallback);
                else
                    InterDivVars_fallback = {};
                end
                
                % call the solver setting `tspan = tspan_fallback`
                [t_fallback, y_fallback] = SOLVER_FALLBACK(ode, tspan_fallback, y0_fallback, options, varargin{:}, InterDivVars_fallback{:});
                
                % write the solution: `2:end` because the last time point of
                % previous iteration = first time point of current iteration
                nt0 = ODEWRITE( t_fallback(2:end), y_fallback(2:end,:), wrapopt.res, nt0);
                                
                if abs(t_fallback(end) - tspan_fallback(end)) > eps
                    error('Fallback solver failed too.')
                end
                
                if isa(wrapopt, 'handle')
                    wrapopt_restart = copy(wrapopt);
                else
                    wrapopt_restart = wrapopt;
                end
                
                % now restart the ode wrapper with the original solver
                
                % update tspan and initial condition
                tspan_restart = T{i}(i2:end);
                
                y0_restart = y_fallback(end,:)';
                
                % append solution and checkpoints
                wrapopt_restart.append = true;
                wrapopt.checkpoints_last_t_index = nt0c;
                
                [~, y0, nt0, nt0c] = odeWrapper(wrapopt_restart, SOLVER, ode, tspan_restart, y0_restart, options, varargin{:} );
            else
                
                if ~wrapopt.save_checkpoints || isempty(ind)
                    error(struct( ...
                        'identifier', 'odeWrapper:solverFailed', ...
                        'message', sprintf('The solver failed at time t = %g. No checkpoints have been saved', t(end))));
                else
                    error(struct( ...
                        'identifier', 'odeWrapper:solverFailed', ...
                        'message', sprintf('The solver failed at time t = %g. Last saved checkpoint at t = %g.', t(end), t(ind(end)))));
                end
                
            end
        end
        
        clear t y
        
    end

    % varargout{1} = wrapopt.res;
    varargout{1} = statistics;
    varargout{2} = y0;
    varargout{3} = nt0;
    varargout{4} = nt0c;
        
    
end

    function nt0 = matlab_matfile( T, Y, res, nt0 )
        
        % get size of the solution in this step
        %   nt  = number of time steps
        %   neq = dimensionality of the solution (number of equations)
        [nt, neq] = size(Y);

        % compute the index-wise interval in the common solution variable
        interval = nt0 + (1:nt);
        
        res.t.t(interval,1) = double(T);
        res.y.y(interval,1:neq) = double(Y);
        
        % update the last index of the yet written solution
        nt0 = interval(end);
        
    end
    function nt0 = matlab_fwrite( T, Y, res, ~ )
        
        res.t.fwrite(double(T));
        res.y.fwrite(double(Y));
        
        nt0 = [];
        
    end
    function nt0 = c_posix( T, Y, res, ~ )
        % http://undocumentedmatlab.com/blog/explicit-multi-threading-in-matlab-part3
        
        save_Posix(res.t.path, 'ab', double(T));
        save_Posix(res.y.path, 'ab', double(Y(:)));

        nt0 = [];
    end

end
