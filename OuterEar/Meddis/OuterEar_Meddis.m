function [ExtEarPressure, B, A, gainScalars] = OuterEar_Meddis(inputPressure, sampleRate, args)

arguments
    inputPressure (1,:) double
    sampleRate (1,1) double
    
    % outer ear resonances band pass filter      [gain order lp hp]
    args.externalResonanceFilters (:,4) double = [10 1 1000 4000];
end
    
Nyquist = sampleRate/2;

n_filters = size(args.externalResonanceFilters, 1);

% details of external (outer ear) resonances
gaindBs = args.externalResonanceFilters(:,1);
filterOrder = args.externalResonanceFilters(:,2);
lowerCutOff = args.externalResonanceFilters(:,3);
upperCutOff = args.externalResonanceFilters(:,4);

gainScalars = 10.^(gaindBs/20);

% external resonance coefficients
B = cell(n_filters,1);
A = cell(n_filters,1);

[Signal_Toolbox_licenceAvailable, ~] = license('checkout','Signal_Toolbox');

for j = 1:n_filters
    
    n = filterOrder(j);
    Wn = [lowerCutOff(j) upperCutOff(j)] / Nyquist;
    
    if Signal_Toolbox_licenceAvailable
        [b, a] = butter(n, Wn);
    else
        [b, a] = butter_no_toolbox(n, Wn, 'bandpass');
    end

    B{j} = b;
    A{j} = a;
end

FilterBndry = cell(n_filters,1);

if isempty(inputPressure)
    ExtEarPressure = [];
else
    y = inputPressure;
    % y = zeros(size(inputPressure));

    for j = 1:n_filters
        % any number of resonances can be used
        
        [x, FilterBndry{j}] = filter( ...
                B{j}, A{j}, inputPressure, FilterBndry{j});

        x = x * gainScalars(j);

        % This is a parallel resonance so add it
        y = y + x;
    end

    ExtEarPressure = y; % pressure at tympanic membrane
end

end
