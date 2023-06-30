function [ runopt, setStructure, stim_struct ] = ...
    findResult( ...
        PART, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, params )
%FINDRESULT Creates settings hash and tries to find saved results

arguments
    PART (1,:) char
    
    stimulus
    topt
    midopt
    mechopt
    mnaopt
    antopt
    hhopt
    runopt
    opt
    
    params.delete_on_error (1,1) logical = true
    params.result_base_dir (1,:) char = ''
end

resultFound = false;

if isempty(params.result_base_dir)
    RES_BASE = opt.resDir;
else
    RES_BASE = params.result_base_dir;
end

%% SWITCH PARAMETERS

required_files_warn = {};
subfolder = '';

switch PART
     % === MID =============================================================
    case 'middleear' 

        RES_SUB = 'mid';
        
        required_files = {'SUCCESS', 't.mat', 'y.mat'};
        
        required_files_warn = {'options.mat', 'settings.mat'};
    
        required_files_OR = {};

    case 'me_mna_dae_ic'

        RES_SUB = PART;

        required_files = {'SUCCESS', 'dae_ic.mat'};
        required_files_warn = {};
        required_files_OR = {};
    
    % === MECH ============================================================
    case 'mechanical'

        RES_SUB = 'mech';
        
        switch mechopt.identifier
            case 'Vetesnik'
                required_files_analysis = { ...
                    fullfile('analysis', 'BMx.mat'), ...
                    fullfile('analysis', 'BMv.mat'), ...
                    fullfile('analysis', 'TMx.mat'), ...
                    fullfile('analysis', 'TMv.mat')};
                switch mechopt.save_method
                    case 'matlab_matfile'
                        required_files = {'SUCCESS', 'statistics.mat', 'mech.mat', 'BM.mat', 'settings.mat', 'cp_t.mat', 'cp_y.mat'};
                    case 'matlab_fwrite'
                        required_files = {'SUCCESS', 'statistics.mat', 'mech.mat', 'BM.mat', 'settings.mat', 'cp_t.mat', 'cp_y.mat'};
                    case 'c_posix'
                        required_files = {'SUCCESS', 'statistics.mat', 'mech.mat', 'BM.mat', 'settings.mat', 'cp_t.mat', 'cp_y.mat'};
                    otherwise
                        error('Save method %s is not recognised option.', mechopt.save_method)
                end
            case 'Verhulst'
                required_files = {'SUCCESS', 'BM.mat'};
        end
        
        required_files_OR = {};
        
        
    % === MNA =============================================================
    case 'oc_mna_circuit_ic'

        RES_SUB = PART;

        required_files = {'SUCCESS', 'circuit_ic.mat'};
        required_files_warn = {};
        required_files_OR = {};
        
    case 'oc_mna_dae_ic'

        RES_SUB = PART;

        required_files = {'SUCCESS', 'dae_ic.mat'};
        required_files_warn = {};
        required_files_OR = {};
    
    case 'oc_mna' 

        RES_SUB = 'mna';
        
        required_files = {'SUCCESS'};
        
        required_files_warn = {'options.mat', 'settings.mat'};
        
        required_files_OR = { ...
            't.mat', 'Time.mat'; ...
            'Y.mat', 'Volt.mat'; ...
            'Y.mat', 'Curr.mat'};

    case 'oc_mna_analysis' 

        RES_SUB = 'mna';
        
        subfolder = 'analysis';
        
        required_files = {};
        required_files_OR = {};
        
        
    % === SYNAPSE =========================================================
    case 'synapse_ic'

        RES_SUB = PART;

        required_files = {'SUCCESS', 'synapse_ic.mat'};
        required_files_warn = {};
        required_files_OR = {};

    case 'synapse' 

        RES_SUB = 'synapse';
        
        % required_files = {'SUCCESS', 'settings.mat'};
        required_files = {'settings.mat'};
        
        required_files_OR = {};
        
    % === NERVE ===========================================================
    case 'nerve'

        RES_SUB = 'nerve';
        
        % required_files = {'SUCCESS', 'settings.mat'};
        required_files = {'settings.mat'};
        
        required_files_OR = {};
        
    % =====================================================================
    otherwise
        error('PART %s is not recognised option.', PART);
end


%% GET RESULT PATH

[ setStr, setHash, setStructure ] = settingsHashAndStr( PART, ...
        [], midopt, mechopt, mnaopt, antopt, hhopt, runopt );

folderName = setHash;
if ~isempty(setStr)
   folderName = [setStr, '_', folderName];
end

resultTopDir = fullfile( RES_BASE, RES_SUB, folderName);

if isempty(stimulus)
    resultPath = fullfile(resultTopDir, subfolder);
else
    [stim_name, stim_hash, stim_struct ] = settingsHashAndStr( ...
            'stimulus', stimulus);

    stim_folder_name = [stim_name, '_', stim_hash];
    resultPath = fullfile( resultTopDir, stim_folder_name, subfolder);
    
    runopt.hash.stimulus = stim_hash;
    runopt.name.stimulus = stim_name;
end

%% CHECK
% * check if result dir exist
% * check if required files are present

% * handle time

%   => return TRUE
% * otherwise
%   => delete files, return FALSE

if strcmp(runopt.purge, PART)
    warning('runopt.purge = true => deleting the folder content ... \n %s', resultPath);
    removeDirAndContents(resultPath);
    error('Forced stop after purging a folder')
    
elseif exist(resultPath, 'dir') && dirsize(resultPath) > 0

    if exist(fullfile(resultPath, 'FAILURE'), 'file')
        warning('FAILURE file found in %s: \n%s', fullfile(resultPath, 'FAILURE'), fileread(fullfile(resultPath, 'FAILURE')))
    end
    
    files_found = zeros(size(required_files));
    if ~isempty(required_files)
        required_files_paths{length(required_files)} = [];
    end

    for j = 1:length(required_files)
        required_files_paths{j} = fullfile(resultPath, required_files{j});

        files_found(j) = exist(required_files_paths{j}, 'file');

    end
    
    for j = 1:length(required_files_warn)
        p = fullfile(resultPath, required_files_warn{j});
        if ~exist(p, 'file')
            warning('Expected file %s, does not exist. Still, we will try to continue.', p);
        end

    end
    
    if ~isempty(required_files_OR)
        files_found_OR = zeros(size(required_files_OR));
%         required_files_paths_OR{[size(required_files_OR)]} = [];

        for i = 1:size(files_found_OR, 1)
            for j = 1:size(files_found_OR, 2)            
                required_files_paths_OR{i,j} = fullfile(resultPath, required_files_OR{i,j});
                files_found_OR(i,j) = exist(required_files_paths_OR{i,j}, 'file');
            end
        end

        files_found_OR_processed = any(files_found_OR, 2);
    else
        files_found_OR_processed = [];
    end

    if any([files_found(:); files_found_OR_processed(:)] == 0)
        
        missingFiles = required_files_paths(~logical(files_found));
        %files_to_delete = required_files_paths(logical(files_found));
        
        s1 = sprintf(' * folder ''%s'' FOUND\n', resultPath);
        s2 = sprintf(' * file ''%s'' MISSING\n', missingFiles{:});                
        
        if params.delete_on_error == true
            warning('Result ''%s'' loading problem! \n %s %s Deleting the folder content ... \n', PART, s1, s2);
        
            removeDirAndContents(resultPath);
        end
    else
        resultFound = true;
    end
end

%% UPDATE RUNOPT

switch PART
    case 'middleear'
        name = 'mid';
    case 'mechanical'
        name = 'mech';
    case 'electrical'
        name = 'mna';
    otherwise
        name = PART;
end

runopt.found.(name) = resultFound;
runopt.path.(name) = resultPath;
runopt.hash.(name) = setHash;

% runopt.path.([name, '_common']) = fullfile(resultTopDir, 'common');
% 
% if ~exist(runopt.path.([name, '_common']), 'dir')
%     mkdir(runopt.path.([name, '_common']))
% end


%% 
    function removeDirAndContents(path)
%         try
            if exist(path,'dir')
                [status, msg] = rmdir(path, 's'); % 's' => remove all contents
                if status ~= 0
                    error('Folder %s not removed. System error %d, %s', path, status, msg)
                end
            end
%         catch e
%             [st, msg] = cmd_rmdir(path);
%             rethrow(e)
%         end
    end

end
