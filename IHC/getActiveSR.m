function [ SR, nSR ] = getActiveSR( antopt )
%GETACTIVESR 
SR = [];


if isa(antopt.fiber, 'cell')
    
    SR = antopt.fiber;
      
elseif isa(antopt.fiber, 'struct')

    allSRs = fields(antopt.fiber);

    for i = 1:length(allSRs)

        if antopt.fiber.(allSRs{i})
            SR = [SR, allSRs(i)];
        end
    end
end

nSR = length(SR);

end

