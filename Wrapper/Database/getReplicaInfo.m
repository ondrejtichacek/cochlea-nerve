function [ replicas ] = getReplicaInfo( PATH, NAME, modelname )
%GETREPLICAINFO
%
%   PATH ... folder to check
%   NAME ... name of file with placeholder signs (eg. *)
%

empty_replica = struct( ...
        'dir', [], ...
        'data', [], ...
        'exists', [], ...
        'tspan', [], ...
        'potential_action', [], ...
        'action', [], ...
        'seed', [], ...
        'rng', []);

switch modelname
    case 'synapse'
        filename = 'synapseResult.mat';
    case 'nerve'
        filename = 'nerveResult.mat';
    otherwise
        error('unexpected option %s', modelname);
end

if ischar(NAME)
    files = dir(fullfile(PATH,NAME));
    names = {files.name};
elseif iscell(NAME)
    names = NAME;
    for i = 1:length(NAME)
        if ~exist(fullfile(PATH, NAME{i}), 'dir')
            names{i} = {};
        end
    end
end

if ~isempty(names)
    replicas = empty_replica;
    
    for i = 1:length(names)
        
        if isempty(names{i})
            replicas(i).dir = fullfile(PATH, NAME{i});
            replicas(i).exists = false;
        else
            replicas(i).dir = fullfile(PATH, names{i});
            if exist(fullfile(replicas(i).dir, filename), 'file')
                replicas(i).exists = true;
                res = load(fullfile(replicas(i).dir, filename));
                switch modelname
                    case 'synapse'
                        res.synapseResult.PATH = replicas(i).dir;
                        replicas(i).data = res.synapseResult;
                    case 'nerve'
                        res.nerveResult.PATH = replicas(i).dir;
                        replicas(i).data = res.nerveResult;
                    otherwise
                        error('Unrecognised model name %s', modelname)
                end
                replicas(i).tspan = tspanOpt( ...
                    't0', Time( ...
                        replicas(i).data.ft.t(1,1), ...
                        replicas(i).data.ft.unit), ...
                    'tf', Time( ...
                        replicas(i).data.ft.t(end,1), ...
                        replicas(i).data.ft.unit));
                replicas(i).seed = replicas(i).data.seed; 
                replicas(i).rng = replicas(i).data.rng_state;
            else
                replicas(i).exists = false;
            end
        end
    end
    
    if replicas(i).exists == false && isdir(replicas(i).dir)
        srmdir(replicas(i).dir, 's')
    end
    
else
    replicas = [];
end
end
