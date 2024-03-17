function [ out, stat_time ] = ReplicationAverage( FILE_HANDLE, ...
    time, PEAK_PROPERTIES, runopt, antopt, hhopt, memopt, SignalPhase, ...
    Nbins, analysis_span_ind, variables, variable_names, args )
%REPLICATIONAVERAGE 
arguments
    FILE_HANDLE
    time
    PEAK_PROPERTIES
    runopt
    antopt
    hhopt
    memopt
    SignalPhase
    Nbins
    analysis_span_ind    
    variables
    variable_names

    args.do_find_peaks (1,:) logical = true
end
%% INIT

% time = time.ms;

S = FILE_HANDLE{1};
time = S.t.s;
time_orig = time;

% find first time index of analysis
% analysisStartTime = runopt.analysisStart.ms;
% analysisStartTimeIndex = min([length(time), find(analysisStartTime <= time,1)]);

% restrict variables to analysis time
time = time(analysis_span_ind(1):analysis_span_ind(2));
SignalPhase = SignalPhase(analysis_span_ind(1):analysis_span_ind(2));

[Signal_Toolbox_licenceAvailable, ~] = license('checkout','Signal_Toolbox');

%% MEMORY REQIUIREMENTS

% L = antopt.numberOfRepetitions * antopt.nSR;

% L = numel(FILE_HANDLE);

% mLim = memoryLimit( memopt );
% [mReq, varProp] = memoryRequirement( FILE_HANDLE, VAR_NAME );

% M = max(1, floor(max(varProp.size(:))*mLim/mReq));

% r = 0:M:L;
% if r(end) ~= L
%     r(end+1) = L;
% end

%% INIT ONLINE CALCULATION OF MEAN & VAR, MAX & MIN

if Signal_Toolbox_licenceAvailable == 1
    PeakFields = PeakStatistics([], [], []);
end

%% CALCULATE MEAN, VAR, MAX, MIN, ...

k = zeros(size(variables));

% p = gcp();
t = tic;
for j = 1:numel(FILE_HANDLE)
    
    S = FILE_HANDLE{j};

    t_preload = tic;
    S = S.preload();
    dbprintf('preload took %s\n', disp_toc(toc(t_preload)))
    
    assert(all(time_orig == S.t.s), 'Different time vectors across files')

    sr = S.SR;

%     X{length(variables)} = [];
%     parfor v = 1:length(variables)
%         X{v} = S.(variables{v});
%     end
    
    for v = 1:length(variables)

        vname = variable_names{v};

        % Since S is a results object, calling S.(variables{v}) will
        % load the variable from the file or compute it on the fly.
        % Therefore, we can preload (precompute) the variable to save
        % time.
        
        N = S.n_dyn_rep; % this line takes a lot of time
        X = S.(variables{v});

        if ~isa(X, 'cell')
            N = 1;
        end

        for ii = 1:N
    
            if isa(X, 'cell')
                x = X{ii};
            else
                x = X;
            end
    
            if isa(x, 'Frequency')
                x = x.Hz;
            elseif ~isa(x, 'double')
                x = double(x);
                warning('Please fix me')
            end
            % x = X{v};
            
            % restrict x to analysis time
            x = x(analysis_span_ind(1):analysis_span_ind(2),:);
            
            n = size(x,2);
            for i = 1:n
                
                if k(v) == 0
                    k(v) = 1;
                    out.(vname)(k(v)) = init_statistics(i, 1, sr, PeakFields, time, Signal_Toolbox_licenceAvailable);
                else
                    % SELECT id FROM ... WHERE v.sr == sr AND v.xind == i
                    kk = find(strcmp({out.(vname).sr}, sr) &  cellfun(@(y) isequal(y, i), {out.(vname).xind}));
                    if numel(kk) > 1
                        error('Max. one result expected');
                    end
                    if isempty(kk)
                        k(v) = numel(out.(vname)) + 1;
                        out.(vname)(k(v)) = init_statistics(i, 1, sr, PeakFields, time, Signal_Toolbox_licenceAvailable);
                    else
                        k(v) = kk;
                    end
                end
                
                out.(vname)(k(v)) = ReplicationStatistics( ...
                    x(:,i), time, out.(vname)(k(v)), PEAK_PROPERTIES, Signal_Toolbox_licenceAvailable, ...
                    'do_find_peaks', args.do_find_peaks(v));
    
    %             % Phase-based statistics
    %             if Signal_Toolbox_licenceAvailable
    %         
    %                 [~, PHI, I, V_phase] = PeriodAverage( ...
    %                         SignalPhase, N, x(:,i));
    %         
    %                 out.(vname)(k(v)).PHASE = V_phase;
    %             end
            end
        end
        
    end
    
%     X{length(variables)} = [];
%     parfor v = 1:length(variables)
%         X{v} = S.(variables{v});
%     end
% 
%     for v = 1:length(variables)
% 
%         pf(v) = parfeval(p, @ReplicationStatistics, 1, ...
%                     X{v}, Time, out.(vname).(sr), PEAK_PROPERTIES, Signal_Toolbox_licenceAvailable);
%     end
% 
%     for v = 1:length(variables)
%         [idx, value] = fetchNext(pf);
%         out{idx}.(sr) = value;
%     end

    S.unload();


    done_frac = j / numel(FILE_HANDLE);
    time_expected = (1-done_frac)*toc(t)/(done_frac);

    dbprintf('Processing file %d/%d, %s expected\n', j, numel(FILE_HANDLE), disp_toc(time_expected))

end
disp_toc(toc(t));


%% Phase-based statistics

if false
    if Signal_Toolbox_licenceAvailable
        for v = 1:length(variables)
            vname = variable_names{v};
            for k = 1:numel(out.(vname))
                [ ~, ~, ~, out.(vname)(k).phase ] = PeriodAverage( SignalPhase, Nbins, out.(vname)(k).MEAN );
            end
        end
    end
end

%% Time-based statistics
if false
    for v = 1:length(variables)
        vname = variable_names{v};
        for k = 1:numel(out.(vname))
            for j = 1:numel(stimulus)
                [ out.(vname)(k).TT{j}, ...
                  out.(vname)(k).C{j}, ...
                  out.(vname)(k).J{j}, ...
                  out.(vname)(k).U{j}, ...
                  out.(vname)(k).I{j} ] = TimeIndexPhaseStat( out.(vname)(k).MEAN, time, Voltage, analysisStartTimeIndex, stimulus(j) );
            end
        end
    end
end
%%

stat_time = time;

end


%% to nahore je zobecneni tohoto
% L = hhopt.numberOfRepetitions * hhopt.nSR;
% 
% % count the number of same SR and  add results
% M = 1000;
% 
% r = 0:M:L;
% if r(end) ~= L
%     r(end+1) = L;
% end
% 
% for j = 1:length(r)-1
%                 
%     S = SynResFile.SynRes(1, (r(j)+1):r(j+1) );
%     
%     
%     for k = 1:length(S)
%     
%         if ~exist('NN', 'var') || ~isfield(NN,S(k).name)
%             NN.(S(k).name) = 1;
%         else
%             NN.(S(k).name) = NN.(S(k).name) + 1;
%         end
% 
%         if ~exist('Cleft', 'var') || ~isfield(Cleft,S(k).name)
%             Cleft.(S(k).name) = S(k).c;
%         else
%             Cleft.(S(k).name) = Cleft.(S(k).name) + S(k).c;
%         end
%         
%         if ~exist('FreePool', 'var') || ~isfield(FreePool,S(k).name)
%             FreePool.(S(k).name) = S(k).q;
%         else
%             FreePool.(S(k).name) = FreePool.(S(k).name) + S(k).q;
%         end
%         
%         if ~exist('ReStore', 'var') || ~isfield(ReStore,S(k).name)
%             ReStore.(S(k).name) = S(k).w;
%         else
%             ReStore.(S(k).name) = ReStore.(S(k).name) + S(k).w;
%         end
%         
%     end
%         
% end
% 
% % get SR names
% SR = fields(NN);
% 
% for i = 1:length(SR)
%     
%     Cleft.(SR{i}) = Cleft.(SR{i}) ./ NN.(SR{i});
%     FreePool.(SR{i}) = FreePool.(SR{i}) ./ NN.(SR{i});
%     ReStore.(SR{i}) = ReStore.(SR{i}) ./ NN.(SR{i});
%     
% end
