function [ status ] = odewaittext( t, ~, flag, tf, ws )
%ODEWAITTEXT
%
% USAGE:
%
% ws = WaitbarStorage();
% ws.num_div = NumDiv;
% 
% OutputFcn = @(t, y, flag) odewaittext(t, y, flag, tspan(end), ws);
% 
% % ~ must correspond with the number of extra arguments in ode15s
% OutputFcnWrapper = @(t, y, flag, ~, ~, ~) OutputFcn(t, y, flag);
% solveropt = odeset(solveropt, 'OutputFcn', OutputFcnWrapper);

if isempty(flag)
    
    printMessageNow = false;
    
    if ( ws.t + ws.print_threshold_progress*tf <= t ) % print message after 5% increase
        printMessageNow = true;
    elseif etime(clock, ws.old_time) >= ws.print_threshold_time  % alternatively print message after 60 sec
        printMessageNow = true;
    end

    if printMessageNow
        
        current_time = clock;
        progress = t(end)/tf;
        progress_increase = (t(end) - ws.t)/tf;
        
        time_diff = etime(current_time,ws.old_time);
        remaining_time = time_diff / progress_increase * (1 - progress);
        
        fprintf('%s : %d %% finished, %s remaining at current speed \n', ...
            datestr(now), round(100*progress), disp_toc(remaining_time));
        
        ws.t = t(end);
        ws.old_time = current_time;
        
    end
    
elseif strcmp(flag, 'init')
    if isempty(ws.division_index)
        ws.division_index = 1;
    end
    if ws.division_index == 1
        ws.t = t(1);
        ws.old_time = clock;
        ws.start_time = clock;
    end
    
elseif strcmp(flag, 'done')
    ws.division_index = ws.division_index + 1;
    if ischar(ws.num_div) || ws.division_index > ws.num_div
        clear ws
    end
    
else
    display(flag)
    error('Unrecognized value of the flag parameter.');
end

status = false;  % to not stop execution

end