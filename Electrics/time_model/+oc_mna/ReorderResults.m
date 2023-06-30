function [ resFiles, statistics ] = ReorderResults( ...
    stimulus, topt, mechopt, mnaopt, runopt, opt, memopt, overwrite )
%REORDERRESULTS Reconstructs the solution variables
arguments
    stimulus
    topt
    mechopt (1,1) mechOpt
    mnaopt (1,1) mnaOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    overwrite (1,1) logical = false
end

%%
Numstacks = mnaopt.Numstacks;

%% Load saved statistics file if exist

f_statistics = fullfile(runopt.path.oc_mna, 'statistics_reconstructResults.mat');

if exist(f_statistics, 'file')
    load(f_statistics, 'statistics');
else
    statistics = [];
end

%%

fmech = fullfile(runopt.path.mech, 'mech.mat');
fBM = fullfile(runopt.path.mech, 'BM.mat');

fMNA = struct( ...
    'fT', fullfile(runopt.path.oc_mna, 'Time.mat'), ...
    'fVolt', fullfile(runopt.path.oc_mna, 'Volt.mat'), ...
    'fCurr', fullfile(runopt.path.oc_mna, 'Curr.mat'), ...
    'fChannels', fullfile(runopt.path.oc_mna, 'Channels.mat'));

if strcmp(mechopt.integration, 'electric')
    fMNA.fMech = fullfile(runopt.path.oc_mna, 'Mech.mat');
end

resFiles.fmech = matfile(fmech);
resFiles.fBM = matfile(fBM);

if overwrite
    delete(fMNA.fT)
    for f = string(fields(fMNA)')
        delete(fMNA.(f))
    end
end

% First check if the required tspan has been already transformed
RequiredTspanTransformed = checkTransformation(topt, fMNA);

fexist = structfun(@(x) exist(x, 'file'), fMNA, 'UniformOutput', true);
if all([RequiredTspanTransformed; ...
        ~overwrite; ...
        fexist])
    
    for f = string(fields(fMNA)')
        resFiles.(f) = matfile(fMNA.(f));
    end
    return

end


%% Transform Results

disp_title(sprintf('transforming MNA results'));

tic

lock_file_cleanup = handleLockFile(runopt.path.oc_mna, 'action', 'create');

% Time from the new file
ft = fullfile(runopt.path.oc_mna, 't.mat');

if ~exist(ft, 'file')
    dir(runopt.path.oc_mna)
    error('Result file %s not found', ft)
end

resFile_t = matfile(ft);
Time = resFile_t.t;
unit_t = resFile_t.unit;
numTime = length(Time);

% Stored time
mfT = matfile(fMNA.fT, 'Writable', true);

try
    l = max(matfileWhosByName( mfT, 'T', 'size'));
    if isempty(l)
        l = 0;
    else
        l = l - 1;
    end
catch ME
    switch ME.identifier
        case 'MATLAB:MatFile:NoFile'
            l = 0;
        otherwise
            rethrow(ME)
    end
end

% Apped new Time to the old
i0 = l + 1;
i1 = l + numTime;

mfT.T(i0:i1,1) = Time;
mfT.unit = unit_t;

% Now handle other variables
fY = fullfile(runopt.path.oc_mna, 'Y.mat');

mReq = memoryRequirement(fY, 'y');
for f = string(fields(fMNA)')
    if exist(fMNA.(f), 'file')
        mReq = mReq + memoryRequirement(fMNA.(f));
    end
end

mLim = memopt.maxVarMem;

processPerPartes = mReq > mLim;

fprintf('memory required for simple load %s, available memory %s\n', ...
    sprintfBytes('%g %s', mReq, 'B'), ...
    sprintfBytes('%g %s', mLim, 'B'));

if processPerPartes

    mReq = mReq/2;
    numTimeToLoad = max(1, min(numTime, floor(numTime*mLim/mReq)));

    timeLoadPoints = unique([0:numTimeToLoad:numTime,numTime]);

    resY = matfile(fY);

    if strcmp(mechopt.integration, 'electric')
        mfMech = matfile(fMNA.fMech,'Writable',true);
    end
    mfCurr = matfile(fMNA.fCurr,'Writable',true);
    mfVolt = matfile(fMNA.fVolt,'Writable',true);
    mfChannels = matfile(fMNA.fChannels,'Writable',true);

    fprintf('loading up to %d time steps at once\n', numTimeToLoad);

else
    timeLoadPoints = [0, numTime];

    if l == 0

        ZERO = zeros(numTime,Numstacks);

        tmp = [cellfun(@(x) sprintf('N%d', x), num2cell(1:mnaopt.circuit.num_node), 'UniformOutput', false); ...
                repmat({ZERO}, 1, mnaopt.circuit.num_node)];
        
        mfVolt = struct(tmp{:});

        mfCurr = struct(...
                    'IHC', ZERO, ...
                    'OHC', ZERO, ...
                    'IHC_MET', ZERO, ...
                    'OHC_MET', ZERO, ...
                    'SV', ZERO ...
                    );
                
        mfChannels = struct(...
                    'IHC', ZERO ...
                    );

        if strcmp(mechopt.integration, 'electric')
            ZERO_MECH = zeros(numTime,mechopt.Numstacks);

            mfMech = struct(...
                        'BMx', ZERO_MECH, ...
                        'TMx', ZERO_MECH ...
                        );
        end

    else
        ZERO = zeros(numTime-1,Numstacks);

        mfCurr = load(fMNA.fCurr);
        mfVolt = load(fMNA.fVolt);
        mfChannels = load(fMNA.fChannels);

        tmp = fields(mfCurr);
        for i = 1:numel(tmp)
            mfCurr.(tmp{i}) = [mfCurr.(tmp{i}); ZERO];
        end

        tmp = fields(mfVolt);
        for i = 1:numel(tmp)
            mfVolt.(tmp{i}) = [mfVolt.(tmp{i}); ZERO];
        end
        
        tmp = fields(mfChannels);
        for i = 1:numel(tmp)
            mfChannels.(tmp{i}) = [mfChannels.(tmp{i}); ZERO];
        end

        if strcmp(mechopt.integration, 'electric')
            ZERO_MECH = zeros(numTime-1,mechopt.Numstacks);

            mfMech = load(fMNA.fMech);

            tmp = fields(mfMech);
            for i = 1:numel(tmp)
                mfMech.(tmp{i}) = [mfMech.(tmp{i}); ZERO_MECH];
            end
        end
    end
end

for i = 2:length(timeLoadPoints)

    ts = timeLoadPoints(i-1) + 1; % starting load point
    te = timeLoadPoints(i); % ending load point

    if processPerPartes
        y = resY.y(ts:te,:);
        ordering = resY.ordering;
        size_info = resY.size_info;
    else
        load(fY, 'y', 'ordering', 'size_info');
    end
    
    % inverse permutations
    Q(ordering.q) = 1:numel(ordering.q);
    y = y(:, Q);

    voltage_names = cellfun(@(x) sprintf('N%d', x), num2cell(1:mnaopt.circuit.num_node), 'UniformOutput', false);
                
    for j = 1:numel(voltage_names)
        name = voltage_names{j};
        mfVolt.(name)(l+(ts:te),1:Numstacks) = getVoltage(y, name, mnaopt);
    end
    
    current_names = {'IHC', 'OHC', 'IHC_MET', 'OHC_MET', 'StV'};
    for j = 1:numel(current_names)
        name = current_names{j};
        mfCurr.(name)(l+(ts:te),1:Numstacks) = getCurrent(y, name, mnaopt);
    end
    
    channels = mnaopt.channels;
    for j = 1:numel(channels)
        mfChannels.(channels(j).name)(l+(ts:te),1:Numstacks) = getChannel(y, channels(j), mnaopt);
    end

    if strcmp(mechopt.integration, 'electric')

        deps_consts.BMx = mechopt.NMkonst_BM;
        deps_consts.BMv = mechopt.NMkonst_BM;
        deps_consts.TMx = mechopt.NMkonst_OHC_cilia;
        deps_consts.TMv = mechopt.NMkonst_OHC_cilia;

        parts = {'BMx', 'BMv', 'TMx', 'TMv'};
        for ii = 1:numel(parts)
            i1 = size_info.(['mech_', parts{ii}]).start;
            i2 = size_info.(['mech_', parts{ii}]).end;

            mfMech.(parts{ii})(l+(ts:te),1:mechopt.Numstacks) = ...
                deps_consts.(parts{ii}) * y(:, i1:i2);
        end
    end
    
    fprintf('%d/%d time points loaded\n',te,numTime);

end

if ~processPerPartes
    if strcmp(mechopt.integration, 'electric')
        save(fMNA.fMech, '-struct', 'mfMech', '-nocompression', '-v7.3');
    end
    save(fMNA.fCurr, '-struct', 'mfCurr', '-nocompression', '-v7.3');
    save(fMNA.fVolt, '-struct', 'mfVolt', '-nocompression', '-v7.3');
    save(fMNA.fChannels, '-struct', 'mfChannels', '-nocompression', '-v7.3');
end

clear y resY mfVolt mfCurr mfChannels

fprintf('transformed results saved to files:\n');
for f = string(fields(fMNA)')
    fprintf('\t%s\n',fMNA.(f));
end

statistics.t_transform = toc;
disp_toc(statistics.t_transform);

for f = string(fields(fMNA)')
    resFiles.(f) = matfile(fMNA.(f));
end

% delete original results to save space
delete(ft);
delete(fY);

%% Saved statistics file

save(f_statistics, 'statistics');

%%

assert(checkTransformation(topt, fMNA), 'Required tspan was not transformed correctly. Did the calculation converge?')

end

function [ RequiredTspanTransformed ] = checkTransformation(topt, fMNA)

    fexist = structfun(@(x) exist(x, 'file'), fMNA, 'UniformOutput', true);

    RequiredTspanTransformed = false;
    if all(fexist)
        mfT = matfile(fMNA.fT);
        T = mfT.T;
        if find(T <= topt.total.t0.(mfT.unit), 1)
            if find(T >= topt.total.tf.(mfT.unit), 1)

                s1 = max(matfileWhosByName(mfT, 'T', 'size'));

                variable_lengths_equal = true;
                break_loop = false;

                for f = string(fields(fMNA)')
                    
                    if strcmp(f, 'fT')
                        continue
                    end

                    if break_loop == true, break, end

                    mf = matfile(fMNA.(f));
                    
                    vars = who(mf);
                    for i = 1:numel(vars)
                        s2 = matfileWhosByName( mf, vars{i}, 'size');

                        if s1 ~= s2(1)
                            variable_lengths_equal = false;
                            break_loop = true;
                            break
                        end
                    end
                end

                if variable_lengths_equal == true
                    RequiredTspanTransformed = true;
                end
            end
        end
    end

end