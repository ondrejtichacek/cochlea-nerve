function [ SynResFile, runopt, antopt ] = handleSynapseResult( Voltage, Time, stimulus, topt, midopt, mechopt, mnaopt, antopt, runopt, paropt, opt, memopt )
%HANDLESYNAPSERESULT Loads stored synapse simulation or executes new

numRepDesired = antopt.numberOfRepetitions;

seeds = zeros(1, numRepDesired, 'int32');

rng('shuffle','twister')

for i = 1:numRepDesired
    while true
        seed = randi(intmax);
        res_dir = fullfile(runopt.path.synapse, sprintf('syn_%x', seed));
        if ~exist(res_dir,'dir') && ~any(seeds == seed)
            seeds(i) = seed;
            break
        end
    end
end


[ SynResFile, runopt, antopt, ~ ] = handleResult( ...
    'synapse', 'any', seeds, numRepDesired, ...
    {Voltage, Time}, {}, ...
    stimulus, topt, midopt, mechopt, mnaopt, antopt, [], runopt, paropt, opt, memopt);


end
