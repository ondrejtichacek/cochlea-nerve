function [ ResFile, runopt, antopt, hhopt ] = handleResult( PART, flag, seeds, numRepDesired, const_src_data, variable_src_data, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, paropt, opt, memopt)
%HANDLESYNAPSERESULT Loads stored synapse/nerve simulation or executes new

% [runopt, setStructure] = findResult( PART, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt );
[runopt, setStructure] = findResultExtended(opt.get_extra_dirs(), PART, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt);

% Setting up the synapse or the nerve function
% for evaluation
switch PART
    case 'synapse'
        P = 'syn';
        fs = antopt.samplingFrequency;
        model = @(rep) setSynapseResult( rep, ...
            const_src_data{:}, topt, antopt, runopt, paropt, memopt, ...
            antopt.fiber, setStructure ...
            );
    case 'nerve'
        P = 'nrv';
        fs = hhopt.samplingFrequency;
        model = @(rep) setNerveResult( rep, ...
            const_src_data{:}, stimulus, topt, mechopt, mnaopt, antopt, hhopt, ...
            opt, runopt, memopt, paropt, setStructure ...
            );
    otherwise
        P = PART;
end

switch flag
    case 'any'
        fnames = sprintf('%s_*', P);
    case 'req'
        assert(numRepDesired == numel(seeds), 'Missmatch in number of replications')
        fnames{numRepDesired} = [];
        for i = 1:numRepDesired
            fnames{i} = sprintf('%s_%x', P, seeds(i));
        end
    otherwise
        error('unrecognised flag value %s', flag)
end


all_replicas = getReplicaInfo( runopt.path.(PART), sprintf('%s_*', P), PART);

replicas = getReplicaInfo( runopt.path.(PART), fnames, PART);

[replicas, numRepDone, numRepExtendable, numRepToDo, numRepToExtend ] = getReplicaActions( replicas, numRepDesired, topt, fs );

if ~runopt.found.(PART)
    fprintf('%s result was not found ...\n', PART)
else
    fprintf('%s result was found,\n', PART)

    fprintf(['specifically:\n', ...
                '   * # done      = %d\n', ...
                '   * # exandable = %d\n', ...
                'with following actions:\n', ...
                '   * # to extend = %d\n', ...
                '   * # to create = %d\n'], ...
                numRepDone, numRepExtendable, numRepToExtend, numRepToDo)

    if runopt.recalculate.(PART)
        fprintf('    BUT we force to recalculate it all ...\n')
        delete_all_replicas(all_replicas);
        replicas = getReplicaInfo( runopt.path.(PART), fnames, PART);
        [replicas, numRepDone, numRepExtendable, numRepToDo, numRepToExtend ] = getReplicaActions( replicas, numRepDesired, topt, fs );

    else
        if numRepToDo > 0 || numRepToExtend > 0
            fprintf('    BUT not all replicas were computed (or long enough), \n')
            if runopt.recalculate.replications
                fprintf('    AND we force to calculate replicas at the same time,\n')
                fprintf('    THEREFORE we force to recalculate all ...\n')
                delete_all_replicas(replicas);
                replicas = getReplicaInfo( runopt.path.(PART), fnames, PART);
                [replicas, numRepDone, numRepExtendable, numRepToDo, numRepToExtend ] = getReplicaActions( replicas, numRepDesired, topt, fs );
            else
                fprintf('    THEREFORE we have to calculate (extend) the missing ...\n')
            end
        elseif numRepToDo == 0 && numRepToExtend == 0
            fprintf('    AND all replicas were already computed ...\n')
        end
    end
end

ResFile{max([numRepDesired, numel(replicas)])} = [];

if numRepToDo < numRepDesired

    if runopt.verbose >= 2
        fprintf('Loading %s result...\n\t%s\n', PART, runopt.path.(PART));
    end

    if runopt.verbose >= 3
        load(fullfile(runopt.path.(PART), 'settings.mat'), ...
            'setStructure');
    display(setStructure)
    end
end

ind = struct();

for action = {'create', 'extend', 'use', 'notuse'}
    ind.(action{:}) = strcmp({replicas.action}, action{:});
end

for i = find([ind.create])
    seed = seeds(i);
    res_dir = fullfile(runopt.path.(PART), sprintf('%s_%x', P, seed));
    if exist(res_dir, 'dir')
        srmdir(res_dir, 's')
    end
    replicas(i).seed = seed;
    replicas(i).dir = res_dir;
    if ~isempty(variable_src_data)
        replicas(i).src_data = variable_src_data{i};
    end
end

% Handle replicas to be loaded
ResFile(ind.use) = {replicas(ind.use).data};


% Handle replicas to be created/extended
if numRepToDo + numRepToExtend > 0

    if runopt.no_create_only_load
        error('LOAD_ONLY mode: result not found!');
    end

    N = model(replicas(ind.create | ind.extend));
    
    ResFile(ind.create | ind.extend) = N;

end

ResFile(ind.notuse) = [];

    function delete_all_replicas(REP)
        if ~isempty(REP)
            for ii = 1:length(REP)
                if ~isempty(REP(ii).dir)
                    if exist(REP(ii).dir, 'dir')
                        srmdir(REP(ii).dir, 's');
                    end
                end
            end
        end
    end


end
