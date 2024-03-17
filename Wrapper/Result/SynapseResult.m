classdef SynapseResult < Result
    %SYNAPSERESULT 
    
    properties
        
        seed
        rng_state
        
        name
        SR
        
        n

        n_dyn_rep

        xgrid
        xgrid_ind
        
        ChannelsOpenSS
        CalciumCurrent
        TransmitterRelease
        
        tropt
        
        protect_PATH_error
        orig_PATH
        
        fy
        ft
        fV
        fy_dyn

        preloaded_files = struct()
        
        size_info
        size_info_dyn
        
        PATH
    end
    
    properties (Dependent)

    end
    
    properties (Transient)
%         fy
%         ft
%         fV
%         
%         PATH

    end
    
    methods (Static)
       function obj = loadobj(s)
           obj = s;
           obj.PATH = [];
           obj.fy = [];
           obj.ft = [];
           obj.fV = [];
           obj.fy_dyn = [];
       end
       
       function i = index_into_variable(varargin)
            if isempty(varargin)
                i = {':',':'};
            else
                i = varargin;
            end
       end
    end
    
    methods
        
        function synapseResult = SynapseResult(varargin)            
            
        end
        
        function val = y(obj, varargin)
            i = index_into_variable(varargin{:});
            
            val = obj.fy.y(i{:});
            
        end
        
        function val = V(obj, varargin)                        
            i = obj.index_into_variable(varargin{:});
            
            val = Voltage(obj.fV.V(i{:}), obj.fV.unit);
            
        end
        
        function val = vesicle_release_events(obj)
            
            val = obj.fy.vesicle_release_events;
        end

        function val = vesicle_release_events_dyn(obj)
            
            val = obj.fy_dyn.vesicle_release_events;
        end
        
        function val = m(obj, varargin)
            
            si = obj.size_info.CaV13_channels_fraction;
            
            % assert(obj.n == 1);
            % ind = 1:obj.n;
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});
            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end
        
        function val = Ca_blocked(obj, varargin)
            
            si = obj.size_info.CaV13_channels_blocked;
            
            % assert(obj.n == 1);
            % ind = 1:obj.n;
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end
        
        function val = I_total(obj, varargin)
            
            si = obj.size_info.Ca_current;
            
            % assert(obj.n == 1);
            % ind = obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            % i{2} = ind(i{2});
            i{2} = ind;

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end

        function val = C(obj, varargin)
            
            si = obj.size_info.Ca_concentration;
            
            % assert(obj.n == 1);
            % ind = obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            % i{2} = ind(i{2});
            i{2} = ind;

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end
        
        function val = q(obj, varargin)
            
            si = obj.size_info.NT_free;
            
            % assert(obj.n == 1);
            % ind = 2*obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            % i{2} = ind(i{2});
            i{2} = ind;

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end

        function val = Cq(obj, varargin)
            
            si = obj.size_info.Ca_concentration;
            
            % assert(obj.n == 1);
            % ind = obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            % i{2} = ind(i{2});
            i{2} = ind;

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});

            si = obj.size_info.NT_free;
            
            % assert(obj.n == 1);
            % ind = 2*obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            % i{2} = ind(i{2});
            i{2} = ind;

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = val .* obj.fy.y(i{:});
        end

        function val = c(obj, varargin)
            
            si = obj.size_info.NT_cleft;
            
            % assert(obj.n == 1);
            % ind = 3*obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end

        function val = w(obj, varargin)
            
            si = obj.size_info.NT_reprocessing;
            
            % assert(obj.n == 1);
            % ind = 4*obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end
        
        function val = c_prot(obj, varargin)
            
            si = obj.size_info.proton_cleft;
            
            % assert(obj.n == 1);
            % ind = 5*obj.n + (1:obj.n);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end

        function val = q_dyn(obj, varargin)
            si = obj.size_info_dyn.NT_free;            
            val = load_dyn_var_2(obj, si, varargin{:});
        end

        function val = c_dyn(obj, varargin)
            si = obj.size_info_dyn.NT_cleft;
            val = load_dyn_var(obj, si, varargin{:});
        end

        function val = w_dyn(obj, varargin)            
            si = obj.size_info_dyn.NT_reprocessing;
            val = load_dyn_var(obj, si, varargin{:});
        end
        
        function val = c_prot_dyn(obj, varargin)
            si = obj.size_info_dyn.proton_cleft;
            val = load_dyn_var(obj, si, varargin{:});
        end

        function val = load_dyn_var(obj, si, varargin)
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            tmp = obj.fy_dyn.y;
            N = numel(tmp);
            val = cell(N,1);
            for j = 1:N
                val{j} = tmp{j}(i{:});
            end

            % this is actually slower
            % val2 = cellfun(@(x) x(i{:}), tmp, 'UniformOutput', false);
        end
        function val = load_dyn_var_2(obj, si, varargin)
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            % i{2} = ind(i{2});
            i{2} = ind;

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            N = numel(obj.fy_dyn.y);
            val = cell(N,1);
            tmp = obj.fy_dyn.y;
            for j = 1:N
                val{j} = tmp{j}(i{:});
            end
        end
        
        function val = CaV13_num_inactivated(obj, varargin)
            
            si = obj.size_info.CaV13_num_inactivated;
            
            % assert(obj.n == 1);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end
        
        function val = CaV13_num_normal(obj, varargin)
            
            si = obj.size_info.CaV13_num_normal;
            
            % assert(obj.n == 1);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end
        
        function val = CaV13_num_burst(obj, varargin)
            
            si = obj.size_info.CaV13_num_burst;
            
            % assert(obj.n == 1);
            
            ind = si.start : si.end;
            
            i = obj.index_into_variable(varargin{:});
            
            i{2} = ind(i{2});

            if numel(obj.n) > 1
                i{3} = 1:numel(obj.n);
            end
            
            val = obj.fy.y(i{:});
        end
        
        function val = m_inf(obj, varargin)
            
            ind = varargin;
            val = obj.ChannelsOpenSS{1}( ...
                    obj.V(ind{:}));
        end
        
        function val = I(obj, varargin)
            
            ind = varargin;
            val = obj.CalciumCurrent( ...
                    obj.V(ind{:}), ...
                    obj.m(ind{:}), ...
                    obj.E_Ca(ind{:}));
        end
        
        function val = k(obj, varargin)
            
            ind = varargin;
            C = obj.C(ind{:});
            if size(C,3) > 1
                val = zeros(size(C));
                for i = 1:size(C,3)
                    val(:,:,i) = obj.TransmitterRelease{i}(C(:,:,i));
                end
            else
                val = obj.TransmitterRelease{1}(C);
            end
        end
        
        function val = E_Ca(obj, varargin)
            ind = varargin;
            
            val = Voltage(Nernst(obj.C(ind{:}), 1.3e-3, 'charge', 2), 'V');
        end
        
        function val = t(obj, varargin)
            
            if isempty(varargin)
                i = {':'};
            else
                i = varargin;
            end
            
            if ndims(obj.ft.t) == 1
                val = Time(obj.ft.t(i{:},:), obj.ft.unit);
            elseif ndims(obj.ft.t) == 2
                val = Time(obj.ft.t(i{:},1), obj.ft.unit);
            end
        end
        
        function PATH = get.PATH(obj)
            if isempty(obj.PATH)
                if obj.protect_PATH_error
                    warning('Please set the `PATH` property after loading the SynapseResult object from file. Usually, it is the folder containing this file.')
                    PATH = obj.orig_PATH;
                    return
                else
                    error('Please set the `PATH` property after loading the SynapseResult object from file. Usually, it is the folder containing this file.')
                end
            end
            PATH = obj.PATH;
        end
        
        function fy = get.fy(obj)
           if isfield(obj.preloaded_files, 'fy')
               fy = obj.preloaded_files.fy;
           else
               if isempty(obj.fy)
                   obj.fy = matfile(fullfile(obj.PATH, 'y.mat'));
               end
               fy = obj.fy;
           end
        end
        
        function ft = get.ft(obj)
           if isfield(obj.preloaded_files, 'ft')
               ft = obj.preloaded_files.ft;
           else
               if isempty(obj.ft)
                   obj.ft = matfile(fullfile(obj.PATH, 't.mat'));
               end
               ft = obj.ft;
           end
        end
        
        function fV = get.fV(obj)
           if isfield(obj.preloaded_files, 'fV')
               fV = obj.preloaded_files.fV;
           else
               if isempty(obj.fV)
                   obj.fV = matfile(fullfile(obj.PATH, 'V.mat'));
               end
               fV = obj.fV;
           end
        end

        function fy_dyn = get.fy_dyn(obj)
           if isfield(obj.preloaded_files, 'fy_dyn')
               fy_dyn = obj.preloaded_files.fy_dyn;
           else
               if isempty(obj.fy_dyn)
                   obj.fy_dyn = matfile(fullfile(obj.PATH, 'y_dyn.mat'));
               end
               fy_dyn = obj.fy_dyn;
           end
        end

        function obj = preload(obj)
            files = {'fy', 'ft', 'fV', 'fy_dyn'};
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

