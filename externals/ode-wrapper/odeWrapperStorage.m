classdef odeWrapperStorage < handle
    
    properties
        
        t           % temporary storage of the variable t
        y           % temporary storage of the variable y
        
        last_checkpoint = 0    % last written checkpoint time
        
        save_checkpoints
        checkpoint_period
        
        nt = 0      % number of time samples in memory
        neq         % number of equations (dimension of y)
        
        maxt        % maximum number of time samples allowed in memory
        
        ft          % t file handle
        fy          % y file handle
        
        cp_t        % t checkpoint file handle
        cp_y        % y checkpoint file handle
        
        output_functions = {};
        
        statistics = struct( ...
            'save_time', 0, ...
            'transfer_time', 0)
            
    end
    
    methods
        function obj = odeWrapperStorage(tspan, y0, wrapopt)
            
            assert(isvector(y0), 'Initial condition must be a vector.')
            
            obj.neq = numel(y0);
            
            % handle maxt
            if ischar(wrapopt.maxt) && strcmp(wrapopt.maxt, 'auto')
                memReq = (obj.neq + 1) * sizeof(class(y0));
                memLim = wrapopt.memopt.maxVarMem;

                wrapopt.maxt = floor(memLim/memReq);
                
                dbprintf('Using %g MB of memory for integration,\ncorresponding to %d time steps.\n', memLim/1024/1024, wrapopt.maxt)
            end
            if numel(tspan) > 2
                obj.maxt = min(numel(tspan), wrapopt.maxt);
            else
                obj.maxt = wrapopt.maxt;
            end
            
            % prealocate arrays
            obj.t = NaN(obj.maxt, 1);
            obj.y = NaN(obj.maxt, obj.neq);
            
            % store file handles
            obj.ft = wrapopt.res.t;
            obj.fy = wrapopt.res.y;
            
            % Handle checkpoints
            if wrapopt.save_checkpoints
                obj.save_checkpoints = wrapopt.save_checkpoints;
                obj.checkpoint_period = wrapopt.checkpoint_period;
                
                obj.cp_t = wrapopt.checkpoints.t;
                obj.cp_y = wrapopt.checkpoints.y;
            end
            
            obj.output_functions{end+1} = @obj.odeSave;
            
        end
        % ----------------------------------------------------------------------
        function status = evalOutputFunctions(obj, t, y, flag)
            % Evaluates all functions stored in obj.output_functions.
            
            status = false;
            for i = 1:numel(obj.output_functions)
                s = obj.output_functions{i}(t, y, flag);
                if s == true
                    warning('Output function %d returned status = true.', i)
                end
                status = status || s;
            end
        end
        % ----------------------------------------------------------------------
        function status = odeSave(obj, t, y, flag)
            
            if isempty(flag)
        
                obj.newIteration(t, y)

            elseif strcmp(flag, 'init')

                obj.newIteration(t(1), y)

            elseif strcmp(flag, 'done')

                obj.finalize()

            else
                display(flag)
                error('Unrecognized value of the flag parameter.');
            end

            status = false;  % to not stop execution
            
        end
        % ----------------------------------------------------------------------
        function newIteration(obj, t, y)
            
            timer = tic;
            
            di = numel(t);
            
            interval = (obj.nt + 1) : (obj.nt + di);
            
            obj.t(interval,1) = t;
            obj.y(interval,:) = y';
            
            obj.nt = obj.nt + di;
            
            obj.statistics.transfer_time = obj.statistics.transfer_time + toc(timer);
            
            if obj.nt >= obj.maxt
                dbprintf(sprintf('Writing to file at time %g', t(end)))
                obj.writeResultToFile();
                obj.nt = 0;
            end
            
            if obj.save_checkpoints
                if abs(obj.last_checkpoint - t(end)) >= obj.checkpoint_period
                    obj.writeCheckpointToFile();
                    obj.last_checkpoint = t(end);
                end
            end
            
        end
        % ----------------------------------------------------------------------
        function finalize(obj)
            
            obj.writeResultToFile();
            
            if obj.save_checkpoints
                obj.writeCheckpointToFile();
            end
            
        end
        % ----------------------------------------------------------------------
        function writeResultToFile(obj)
            
            timer = tic;
            
            if obj.nt > 0
                obj.writeToFile('ft', 'fy', 1:obj.nt, 1:obj.neq);
            end
            
            obj.statistics.save_time = obj.statistics.save_time + toc(timer);
        end
        % ----------------------------------------------------------------------
        function writeCheckpointToFile(obj)
            
            timer = tic;
            
            if obj.nt > 0
                obj.writeToFile('cp_t', 'cp_y', obj.nt, 1:obj.neq);
            end
            
            obj.statistics.save_time = obj.statistics.save_time + toc(timer);
        end
        % ----------------------------------------------------------------------
        function writeToFile(obj, FT, FY, interval, yinterval)
            
            save_to_matfile_nocompression(obj.(FT), 't', obj.t(interval))
            save_to_matfile_nocompression(obj.(FY), 'y', obj.y(interval, yinterval))
            
%             t_saved = false;
%             y_saved = false;
%             
%             ft_exists = exist(obj.(FT).Properties.Source, 'file');
%             
%             if ~ft_exists || isempty(who(obj.(FT), 't'))
%                 
%                 if ft_exists
%                     saveflags = {'-append'}; % v7.3 not included to prevent warning
%                 else
%                     saveflags = {'-v7.3'};
%                 end
%                 
%                 t = obj.t(interval);
%                 
%                 save(obj.(FT).Properties.Source, ...
%                     '-nocompression', saveflags{:}, 't');
%                 
%                 % obj.(FT) = matfile(obj.(FT).Properties.Source, ...
%                 %     'Writable', true);
%                 
%                 t_saved = true;
%             end
%             
%             fy_exists = exist(obj.(FY).Properties.Source, 'file');
%             
%             if ~fy_exists || isempty(who(obj.(FY), 'y'))
%                 
%                 if fy_exists
%                     saveflags = {'-append'}; % v7.3 not included to prevent warning
%                 else
%                     saveflags = {'-v7.3'};
%                 end
%                 
%                 y = obj.y(interval, yinterval);
%                 
%                 save(obj.(FY).Properties.Source, ...
%                     '-nocompression', saveflags{:}, 'y');
%                 
%                 % obj.(FY) = matfile(obj.(FY).Properties.Source, ...
%                 %    'Writable', true);
%                 
%                 y_saved = true;
%             end
%             
%             assert(t_saved == y_saved)
%             
%             if ~t_saved
%                 
%                 nt0 = max(matfileWhosByName(obj.(FT), 't', 'size'));
%                 
%                 assert(nt0 > 0)
% 
%                 int = 1:numel(interval);
% 
%                 obj.(FT).t(nt0 + int, 1) = obj.t(interval);
%                 obj.(FY).y(nt0 + int, yinterval) = obj.y(interval, yinterval);
%             end
        end
    end
end