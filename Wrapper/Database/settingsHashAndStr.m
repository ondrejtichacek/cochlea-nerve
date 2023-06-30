function [ STR, HASH, setStructure ] = settingsHashAndStr( MODEL_PART, stimulus, midopt, mechopt, mnaopt, antopt, hhopt, runopt )
%SETTINGSHASHANDSTR Summary
arguments
    MODEL_PART (1,:) char
    stimulus = []
    midopt = []
    mechopt = []
    mnaopt = []
    antopt = []
    hhopt = []
    runopt = []
end

if ~strcmp(MODEL_PART, ...
        {'stimulus', ...
        'middleear', ...
        'me_mna_dae_ic', ...
        'mechanical', ...
        'oc_mna', ...
        'oc_mna_circuit_ic', ...
        'oc_mna_dae_ic', ...
        'oc_mna_analysis', ...
        'synapse_ic', ...
        'synapse', ...
        'nerve'})
    
    error('Unrecognized model part %s', MODEL_PART);
end

%% Init settings structure and settings string

STR = '';

if strcmp(MODEL_PART, 'stimulus')
       
    switch class(stimulus)
        case 'ZeroSignal'
            stimulus_name = sprintf('zero');
        otherwise
            stimulus_name = stimulus.name;
    end
    
    STR = [STR, stimulus_name]; %, sprintf('%dc', mechopt.Numstacks)];
    
    setStructure.stimulus = [ stimulus.DefiningData ];
    
else
    setStructure = struct( ...
        'CodeVersion', runopt.CodeVersion ...
    );
end

switch MODEL_PART
    case 'stimulus'
        
    case 'middleear'
        % STR = [STR, ''];
        setStructure.middlear = midopt.DefiningData;
       
    case {'me_mna_dae_ic'}
        setStructure.middlear = midopt.DefiningData;
        
    case 'mechanical'
        % STR = [STR, ''];
        setStructure.middlear = midopt.DefiningData;      
        setStructure.mechanical = mechopt.DefiningData;

    case {'oc_mna', 'oc_mna_analysis'}
        % STR = [STR, ''];
        setStructure.middlear = midopt.DefiningData;
        setStructure.mechanical = mechopt.DefiningData;
        setStructure.electrical = mnaopt.DefiningData;

    case {'oc_mna_circuit_ic', 'oc_mna_dae_ic'}
        % STR = [STR, ''];
        tmp = mechopt.DefiningData;
        
        fn = { ...
            'Numstacks', ...
            'identifier', ...
            'specifier', ...
            'integration'};
        
        for i = 1:numel(fn)
            setStructure.mechanical.(fn{i}) = tmp.(fn{i});
        end
        
        setStructure.electrical = mnaopt.DefiningData;
        
        setStructure.electrical = rmfield(setStructure.electrical, 'samplingFrequency');
    
    case {'synapse_ic'}
        % STR = [STR, ''];
        tmp = mechopt.DefiningData;
        
        fn = { ...
            'Numstacks', ...
            'identifier', ...
            'specifier', ...
            'integration'};
        
        for i = 1:numel(fn)
            setStructure.mechanical.(fn{i}) = tmp.(fn{i});
        end
        
        setStructure.electrical = mnaopt.DefiningData;
        
        setStructure.electrical = rmfield(setStructure.electrical, 'samplingFrequency');

        setStructure.synapse = antopt.DefiningData;

        setStructure.synapse.positionsToExcite = antopt.positionsToExcite;

    case 'synapse'
        STR = [STR, sprintf('%s', ...
            antopt.ant)];

        setStructure.middlear = midopt.DefiningData;
        setStructure.mechanical = mechopt.DefiningData;
        setStructure.electrical = mnaopt.DefiningData;
        setStructure.synapse = antopt.DefiningData;

    case 'nerve'
        STR = [STR, sprintf('%s', ...
            hhopt.method)];

        setStructure.middlear = midopt.DefiningData;
        setStructure.mechanical = mechopt.DefiningData;
        setStructure.electrical = mnaopt.DefiningData;
        setStructure.synapse = antopt.DefiningData;
        setStructure.nerve = hhopt.DefiningData;
        
    otherwise
        error('Unrecognised part %s', MODEL_PART)
end


%%

OPT.Method = 'SHA-1';
OPT.Format = 'base32';
HASH = DataHash(setStructure, OPT);

end
