function [ replicas, Done, Extendable, ToDo, ToExtend ] = getReplicaActions( replicas, numRepDesired, topt, fs )
%GETREPLICAACTIONS
%
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

if isempty(replicas)
    replicas = empty_replica;
end
    

if length(replicas) < numRepDesired
    for i = numRepDesired : -1 : (length(replicas) + 1)
        replicas(i) = empty_replica;
        replicas(i).exists = false;
    end
end
    

numRepDone = 0;
numRepExtendable = 0;

TF = Time(zeros(length(replicas),1));

for i = 1:length(replicas)
    
    if replicas(i).exists
    
        tf_found = replicas(i).tspan.tf;
        tf_required = topt.total.tf;

        TF(i) = tf_found;

        if tf_required <= tf_found || abs(tf_required - tf_found) < 2/fs
            % no need to compute anything
            numRepDone = numRepDone + 1;
            replicas(i).tspan_complement = [];
            replicas(i).potential_action = 'use';
        else
            % need to compute the result from `tf_found` to `tf_required`
            numRepExtendable = numRepExtendable + 1;
            replicas(i).tspan_complement = tspanOpt('t0', tf_found, 'tf', tf_required);
            replicas(i).potential_action = 'extend';
        end

    else
        
        TF(i) = Time(-Inf);
        replicas(i).tspan_complement = copy(topt.total);
        replicas(i).potential_action = 'create';
        
    end
end

numRepToDo = max(0, numRepDesired - numRepDone - numRepExtendable);
numRepToExtend = min(numRepExtendable, numRepDesired - numRepDone);

[~,ind] = sort(TF.ms, 'descend');

for j = 1:length(ind)
    i = ind(j);

    if j <= numRepDesired
        replicas(i).action = replicas(i).potential_action;
    else
        replicas(i).action = 'notuse';
    end
end

Done = numRepDone;
ToDo = numRepToDo;
Extendable = numRepExtendable;
ToExtend = numRepToExtend;


end