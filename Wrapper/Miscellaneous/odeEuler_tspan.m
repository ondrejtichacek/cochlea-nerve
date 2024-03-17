function [ numSteps, numSamples, t ] = odeEuler_tspan(tspan, dt)
%ODEEULER
%

t0 = tspan(1);
tf = tspan(2);

if length(tspan) == 2
    
    numSteps = ceil((tspan(2) - tspan(1)) / dt);
    numSamples = numSteps + 1;
    
    t = linspace(t0, tf, numSamples);
    t = t(:);
    
else
    
    t = tspan(:);
    numSamples = length(t);
    numSteps = numSamples - 1;
    
    meanTimeStepDifference = mean(abs(diff(t) - dt));
    
    if meanTimeStepDifference >= dt*1e-3
        warning('Mean difference in actual time step and requested time step is %g\n', meanTimeStepDifference)
    end
    
end

end