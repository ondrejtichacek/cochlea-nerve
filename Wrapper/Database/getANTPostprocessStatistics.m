function [ ANTStatistics ] = getANTPostprocessStatistics( ...
        ANT_fun, ...
        stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt )
%GETANTPOSTPROCESSSTATISTICS

switch ANT_fun
    case 'ANT'
        if max([stimulus.frequency]) > 0
            mechopt.samplingFrequency = 64 * max([stimulus.frequency]);
            mnaopt.samplingFrequency = 64 * max([stimulus.frequency]);
        end
        FUN = @ANT;

    case 'ANT_single'
        FUN = @ANT_single;

    case 'ANT_clamp'
        FUN = @ANT_clamp;

    otherwise 
        error('asd')
end

% [ runopt ] = findResult( 'synapse', stimulus, topt, midopt, mechopt, mnaopt, antopt, [], runopt, opt );
[runopt, ~] = findResultExtended(opt.get_extra_dirs(), 'synapse', stimulus, topt, midopt, mechopt, mnaopt, antopt, [], runopt, opt);

% postprocess_windows = {'full', 'main', 'onset_20ms'};
% postprocess_windows = {'full'};

if isa(stimulus, 'PureTone')
    postprocess_windows = {...
        'full', 'main', ...
        'onset_10ms', 'offset_10ms', ...
        'onset_20ms', 'offset_20ms', ...
        'ss_zero_last_10ms', ...
        'ss_zero', 'ss_end_50ms'};
else
    postprocess_windows = {'full'};
end

doANT = false;

for i = 1:numel(postprocess_windows)
        w = postprocess_windows{i};
        StatFileName = fullfile(runopt.path.synapse, ...
            sprintf('ANTStatistics_%s.mat', w));
    
    if runopt.recalculate.ant_postprocess
        doANT = true;
    end
    
    if ~exist(StatFileName,'file')
        doANT = true;
    end
    
    if exist(StatFileName,'file')
        load(StatFileName, 'num_replicas_synapse');
        if ~exist('num_replicas_synapse', 'var')
            warning('Variable num_replicas_synapse not found in ANTStatistics.mat, recalculating');
            doANT = true;
        else
            if num_replicas_synapse < antopt.numberOfRepetitions
                doANT = true;
            end
        end
    end
end

if doANT
    fprintf('Recalculating ANTPostprocess statistics\n');
    ANTStatistics = FUN( stimulus, topt, midopt, mechopt, mnaopt, antopt, hhopt, runopt, opt, memopt, paropt, [] );
else
    fprintf('Loading ANTPostprocess statistics\n');
    for i = 1:numel(postprocess_windows)
        w = postprocess_windows{i};
        StatFileName = fullfile(runopt.path.synapse, ...
            sprintf('ANTStatistics_%s.mat', w));
        ANTStatistics.(w) = matfile(StatFileName, 'Writable', true);
    end
end

end
