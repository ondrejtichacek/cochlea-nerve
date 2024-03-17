function [ t, y ] = odeEuler( ode, tspan, y0, options, varargin )
%ODEEULER
%

dt = options.TimeStep;
UseOutputFcn = options.UseOutputFcn;
OutputFcn = options.OutputFcn;
OutputFcnEvalInterval = options.OutputFcnEvalInterval;

t0 = tspan(1);
tf = tspan(2);

[numSteps, numSamples, t] = odeEuler_tspan(tspan, dt);

n = size(y0,1);

y = zeros(numSamples,n);

y(1,:) = y0;

if isempty(options.Precompute)
    precomputed_args = [];
else
    precomputed_args = options.Precompute(t);
    varargin{end+1} = [];
end

OutputFcnEvalTime = t0;
if UseOutputFcn
    OutputFcn(t0,y0,'init');
end

for i = 1:numSteps
    
    if ~isempty(precomputed_args)
        varargin{end} = precomputed_args(i,:);
    end

    y(i+1,:) = y(i,:) + dt * ode(t(i+1), y(i,:)', varargin{:})';
    
    if UseOutputFcn
        if t(i+1) - OutputFcnEvalTime >= OutputFcnEvalInterval
            status = OutputFcn(t(i+1),y(i+1,:)',[]);
            if status == true
                break
            end
        end
    end
    
end

if UseOutputFcn
    OutputFcn([],[],'done');
end

end