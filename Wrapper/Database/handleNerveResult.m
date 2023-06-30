function [ NerveResFile, runopt, hhopt ] = handleNerveResult( SynResFile, stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, paropt, opt, memopt )
%HANDLENERVERESULT Loads stored nerve simulation or executes new

numRepDesired = numel(SynResFile);
if numRepDesired ~= antopt.numberOfRepetitions
    error('Missmatch in number of replications')
end

seeds = cellfun(@(x) x.seed, SynResFile);

[ NerveResFile, runopt, ~, hhopt ] = handleResult( ...
    'nerve', 'req', seeds, numRepDesired, ...
    {}, SynResFile, ...
    stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, paropt, opt, memopt);



end
