classdef NerveResult < Result
    %NERVERESULT 
    
    properties
        n

        n_dyn_rep

        Na_max
        seed
        rng_state
        
        xgrid
        xgrid_ind

        name
        SR
        
        calcN_Na
        mh
        mhinit
        
        protect_PATH_error
        orig_PATH

        preloaded_files = struct()
    end
    
    properties (Dependent)

    end
    
    properties (Transient)
        fy
        ft
        fi
        
        PATH

    end
    
    methods (Static)
       function obj = loadobj(s)
           obj = s;
           obj.PATH = [];
           obj.fy = [];
           obj.ft = [];
           obj.fi = [];
       end
    end
    
    methods
        
        function nerveResult = NerveResult(varargin)            
            
        end
        
        function val = y(obj, varargin)                        
            
            if isempty(varargin)
                i = {':',':'};
            else
                i = varargin;
            end            
            
            val = obj.fy.y(i{:});
            
        end
        
        function val = Istim(obj, varargin)                        
            
            if isempty(varargin)
                i = {':',':'};
            else
                i = varargin;
            end            
            
            val = obj.fi.Istim(i{:});
            
        end
        
        function val = V(obj, varargin)
            
            ind = 1:obj.n;
            
            if isempty(varargin)
                i = {':',':'};
            else
                i = varargin;
            end
            
            i{2} = ind(i{2});
            
            N = obj.n_dyn_rep;
            val = cell(N,1);
            tmp = obj.fy.y;
            for j = 1:N
                val{j} = tmp{j}(i{:});
            end
        end
        
        function val = m(obj, varargin)
            
            ind = obj.n + (1:obj.n);
            
            if isempty(varargin)
                i = {':',':'};
            else
                i = varargin;
            end
            
            i{2} = ind(i{2});
            
            if isa(obj.fy.y, 'cell')
                val = cellfun(@(x) x(i{:}), obj.fy.y, 'UniformOutput', false);
            else
                val = obj.fy.y(i{:});
            end

        end
        
        function val = h(obj, varargin)
            
            ind = 2*obj.n + (1:obj.n);
            
            if isempty(varargin)
                i = {':',':'};
            else
                i = varargin;
            end
            
            i{2} = ind(i{2});
            
            if isa(obj.fy.y, 'cell')
                val = cellfun(@(x) x(i{:}), obj.fy.y, 'UniformOutput', false);
            else
                val = obj.fy.y(i{:});
            end
        end
        
        function val = N_Na(obj, varargin)
            
            ind = varargin;
            val = obj.calcN_Na( ...
                    obj.m(ind{:}), ...
                    obj.h(ind{:}), ...
                    obj.Na_max);
            % val = zeros(size(obj.m(ind{:})));
        end
        
        function val = t(obj, varargin)
            
            if isempty(varargin)
                i = {':'};
            else
                i = varargin;
            end
            
            val = Time(obj.ft.t(i{:},:), obj.ft.unit);
        end
        
        function PATH = get.PATH(obj)
            if isempty(obj.PATH)
                if obj.protect_PATH_error
                    warning('Please set the `PATH` property after loading the NerveResult object from file. Usually, it is the folder containing this file.')
                    PATH = obj.orig_PATH;
                    return
                else
                    error('Please set the `PATH` property after loading the NerveResult object from file. Usually, it is the folder containing this file.')
                end
            end
            PATH = obj.PATH;
        end
        
        function fy = get.fy(obj)
           if isfield(obj.preloaded_files, 'fy')
               fy = obj.preloaded_files.fy;
           else
               if isempty(obj.fy)
                   obj.fy = matfile(fullfile(obj.PATH,'y.mat'));
               end
               fy = obj.fy;
           end
        end
        
        function ft = get.ft(obj)
            if isfield(obj.preloaded_files, 'ft')
                ft = obj.preloaded_files.ft;
            else
                if isempty(obj.ft)
                    obj.ft = matfile(fullfile(obj.PATH,'t.mat'));
                end
                ft = obj.ft;
            end
        end
        
        function fi = get.fi(obj)
           if isfield(obj.preloaded_files, 'fi')
               fi = obj.preloaded_files.fi;
           else
               if isempty(obj.fi)
                   obj.fi = matfile(fullfile(obj.PATH,'i.mat'));
               end
               fi = obj.fi;
           end
        end

        function obj = preload(obj)
            files = {'fy', 'ft', 'fi'};
            for i = 1:numel(files)
                f = files{i};
                if isfield(obj.preloaded_files, f)
                    % already preloaded
                else
                    s = obj.(f).Properties.Source;
                    obj.preloaded_files.(f) = load(s);
                end
            end
        end

        function obj = unload(obj)
            obj.preloaded_files = [];
        end
        
%         function v = subsref(this, s)
%             
%             if strcmp(s(1).type, '.') && any(strcmp(s(1).subs, ...
%                     {'V', 'm', 'h'}))
%                 
%                 if numel(s) == 1
%                     q = {':',':'};
%                 elseif numel(s) == 2
%                     q = s(2).subs;
%                 else
%                     error('Multilevel indexing not supported');
%                 end
%                 
%                 switch s(1).subs
%                     case {'V', 'm', 'h'}
%                         switch s(1).subs
%                             case 'V'
%                                 ind = 1:this.n;
%                             case 'm'
%                                 ind = this.n + (1:this.n);
%                             case 'h'
%                                 ind = 2*this.n + (1:this.n);
%                         end
%                                 
%                         if strcmp(q{2},':'), 
%                             q{2} = ind;
%                         else 
%                             q{2} = ind(q{2});
%                         end
%                         
%                         v = this.fy.y(q{:});
%                 end
%                                 
%             else
%                 v = builtin('subsref',this,s);
%             end
%         end
        
        
    end
    
end

