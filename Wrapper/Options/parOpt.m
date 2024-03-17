classdef parOpt < Opt

    properties
        force_num_cores (1,1) double = 0

        useparfor (1,1) logical = true
        useparfeval (1,1) logical = true
        submit_parfeval (1,1) logical = false
        submitAsBatchJob (1,1) logical = false
        
        custom_monitor (1,1) logical = false
    end
    properties (Transient)
        parcluster
        parpool
        useparalleltoolboxinheuristic (1,1) logical = false
        useparalleltoolboxinouterscript (1,1) logical = false
        useparalleltoolbox (1,1) logical = false
    end
    properties (Dependent)
        licenceAvailable
    end
    properties (Dependent, Hidden)
        definingProperties = {};
    end
    properties (Constant, Hidden)
        requiredParameters = {};
    end

    methods
        
        % =================================================================
        % GET FUNCTIONS FOR PARAMETERS
        % -----------------------------------------------------------------
        %
        function lic = get.licenceAvailable(obj)
            lic = true;
            % lic = license('checkout','Distrib_Computing_Toolbox');
        end
        function val = get.useparalleltoolbox(obj)
            if obj.useparalleltoolbox && ~obj.licenceAvailable
                warning('Can''t use parallel toolbox, license is not available');
                val = false;
            else
                val = obj.useparalleltoolbox;
            end
        end
        % -----------------------------------------------------------------
        % =================================================================
        
        function obj = parOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % check if all required properties are set
            obj.checkRequired(opt);
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
            if obj.useparalleltoolbox && obj.submitAsBatchJob
                if obj.useparfor

                    warning('parfor not available when submitted as batch job, setting useparfor = false');        
                    obj.useparfor = false;

                end
                if obj.useparfeval

                    warning('parfeval not available when submitted as batch job, setting useparfeval = false');
                    obj.useparfeval = false;

                end
            end
        end
        function configure_parcluster(obj, num_cores)
            arguments
                obj
                num_cores = []
            end
            [username, computer] = userAndComputerName();
            
            ppname = randstr(8);
            
            if isAurum()
                fn = fullfile('/dev/shm/slurm_jobs', username, getenv('SLURM_JOB_ID'), 'local_cluster_jobs', ppname);
            else
                fn = fullfile('/dev/shm', username, 'local_cluster_jobs', ppname);
            end
            if exist(fn, 'dir')
            else
                mkdir(fn)
            end

            disp(fn)
            
%             if num_cores == 0
%                 profile_name = 'local';
%             else
%                 profile_name = sprintf('local_%d', num_cores);
%             end
            warning('This is workaround now')
            profile_name = 'Processes';

            pc = parcluster(profile_name);
            if pc.NumWorkers == 1
                disp(pc);
                warning('It seems your parcluster profile is not well configured, try to follow https://www.mathworks.com/help/parallel-computing/saveprofile.html and modify the profile');
                % pc = parcluster()
                % pc.NumWorkers = num_cores
                % saveAsProfile(pc, profile_name)
                % ALLPROFILES = parallel.clusterProfiles
            end
            pc.JobStorageLocation = fn;
            
            if ~isempty(num_cores) && num_cores > 0
                assert(pc.NumWorkers == num_cores);
            end
            
            obj.parcluster = pc;            
        end
        function create_cluster_profiles(obj, max_num_cores)
            
            error('Deprecated')
            
            [username, computer] = userAndComputerName();
            all_profiles = parallel.clusterProfiles;
            for num_cores = 1:max_num_cores
                
                profile_name = sprintf('local_%d', num_cores);
                
                if any(strcmp(profile_name, all_profiles))
                    pc = parcluster(profile_name);
                else
                    pc = parcluster();
                end
                
                pc.NumWorkers = num_cores;
                % fn = fullfile('/dev/shm', username, 'local_cluster_jobs');
                % pc.JobStorageLocation = fn;
                
                saveAsProfile(pc, profile_name)
            end
            disp(parallel.clusterProfiles);
        end
        function p = startparpool(obj, num_cores, args)
            arguments
                obj
                num_cores = []
                args.debug (1,1) logical = false
            end

%             p = gcp("nocreate");
%             if isempty(p)
%                 p = parpool('Processes', 64);
%             end
%             warning('This is workaround now')
%             return


            if isempty(num_cores) && obj.force_num_cores
                num_cores = obj.force_num_cores;
            end

            if isempty(obj.parcluster) && obj.force_num_cores
                obj.configure_parcluster(obj.force_num_cores);
            end

            if isempty(obj.parcluster)
                if args.debug == true
                    dbprintf('obj.parcluster is empty, selecting ''Processes'' profile\n');
                end
                pc = 'Processes';
            else
                if args.debug == true
                    dbprintf('obj.parcluster not empty\n');
                    disp(obj.parcluster);
                end
                pc = obj.parcluster;
            end
            p = gcp('nocreate');
            t = tic;
            if isempty(p)
                if isempty(num_cores) || num_cores == 0
                    p = parpool(pc);
                else
                    p = parpool(pc, num_cores);
                end
            else
                if isempty(num_cores)
                    if p.NumWorkers ~= num_cores
                        delete(p)
                    end
                end
            end
            texec = toc(t);
            if args.debug == true
                fprintf('parpool created in %s\n', disp_toc(texec));
            end
            obj.parpool = p;
        end
    end
    
end