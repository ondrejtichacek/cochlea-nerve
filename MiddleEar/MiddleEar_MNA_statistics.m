function [ stats, res_files ] = MiddleEar_MNA_statistics( ...
    stimulus, topt, midopt, runopt, opt, memopt, paropt, do)
%MIDDLEEAR_MNA_STATISTICS
arguments
    stimulus (1,1) stimulusOpt
    topt   (1,1) timeOpt
    midopt (1,1) midOpt
    runopt (1,1) runOpt
    opt    (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt
    
    do.PhaseLagFFT (1,1) logical = true
    do.HilbertPhase (1,1) logical = true
end

% variables to load from saved data
variables = [ ...
    struct( ...
        'id', 'pcoch', ...
        'name', 'cochlea pressure', ...
        'unit', 'Pa')];

assert(numel(variables) == 1); % so far only one variable
    
res_files = ResCacheMiddleEar(stimulus, topt, midopt, runopt, opt, memopt, paropt);

analysis_folder = fullfile(runopt.path.mid, 'analysis');

% Set the analysis files
for i = 1:numel(variables)
    variables(i).analysis_file = fullfile(analysis_folder, sprintf('%s.mat',variables(i).id));
    variables(i).analysis_file_found = exist(variables(i).analysis_file, 'file');
end

all_analysis_files_exist = all([variables.analysis_file_found]);

if runopt.recalculate.mid || ...
        ( ~all_analysis_files_exist || runopt.recalculate.mid_statistics )

    y = res_files.fy.y;
    t = Time(res_files.ft.t, 's');

    nt = numel(t);
    n = 1;

    t1 = find(runopt.analysisStart <= t, 1);
    if isempty(t1)
        error('runopt.analysisStart is larger than the simulation length')
    end
    t2 = length(t);

    for i = 1:numel(variables)
        switch variables(i).id
            case 'pcoch'
                
                pcoch = y(:,midopt.IND);

                REF = pcoch(1, :);
                VAR = pcoch(t1:t2, :);
                
        end
        stats = analyzeSignalME(VAR, REF, t, t1, t2, n, nt, stimulus, do);
        
        if ~exist(analysis_folder, 'dir')
            mkdir(analysis_folder);
        end
        
        save(variables(i).analysis_file, '-struct', 'stats')
    end
else
    for i = 1:numel(variables)
        stats = load(variables(i).analysis_file);
    end
end
end

function stats = analyzeSignalME(VAR, REF, t, t1, t2, n, nt, stimulus, do)

stats.max = max(VAR);
stats.min = min(VAR);

stats.max_rel = max(VAR-REF);
stats.min_rel = min(VAR-REF);

fs = Frequency(200, 'kHz');
% fs = 10*stimulus.frequency;
dt = 1/fs;

T = Time([t(t1).us:dt.us:t(t2).us]', 'us');

if n == 1
    F = griddedInterpolant(t(t1:t2).us, VAR);
    VAR = F(T.us);
else
    F = griddedInterpolant({t(t1:t2).us, 1:n}, VAR);
    VAR = F({T.us, 1:n});
end

assert(size(VAR,1) == length(T), 'Interpolated voltage has wrong shape')

S = stimulus.eval(T);

if do.PhaseLagFFT == true
    x = S;
    y = VAR;

    x = x - mean(x);
    y = y - mean(y);

    L = length(x);
    assert(L == length(y));

    X = fft(x);
    Y = fft(y);        

    f = fs.Hz*(0:(L/2))/L;

    % P2 = abs(X/L);
    % P1 = P2(1:L/2+1);
    % P1(2:end-1) = 2*P1(2:end-1);

    % plot(f,P1) 
    % title('Single-Sided Amplitude Spectrum of S(t)')
    % xlabel('f (Hz)')
    % ylabel('|P1(f)|')

    [~, idx_x] = max(abs(X));
    [~, idx_y] = max(abs(Y));

    % determine the phase difference at the maximum point.
    px = angle(X(idx_x));
    py = angle(Y(idx_y));

    stats.phase_lag_fft = py - px;
    % stats.phase_lag_fft = px - py;

    stats.f_stimulus = f(idx_x);
    stats.f_signal = f(idx_y);

end

if do.HilbertPhase == true
    [Signal_Toolbox_licenceAvailable, ~] = license('checkout','Signal_Toolbox');
    if Signal_Toolbox_licenceAvailable

        stats.phase_lag_hilbert = mean(unwrap(angle(hilbert(VAR))) - unwrap(angle(hilbert(S))));
        % stats.phase_lag_hilbert = mean(unwrap(angle(hilbert(S))) - unwrap(angle(hilbert(VAR))));

    else
        warning('Can''t use hilbert transform without signal processing toolbox. Not including phases in the ouput.')
        stats.phase_lag_hilbert = [];
    end
end
end
