function [ variables, runopt ] = oc_mna_analysis( ...
    resFiles, stimulus, topt, midopt, mechopt, mnaopt, runopt, opt, memopt, paropt, args )
%MNA
arguments
    resFiles
    stimulus 
    topt (1,1) timeOpt
    midopt (1,1) midOpt
    mechopt (1,1) mechOpt
    mnaopt (1,1) mnaOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt
    args.all_analysis_files_must_exist (1,1) logical = false
    args.variables_to_analyze = struct( ...
        'voltage', 'all', ...
        'current', 'all', ...
        'mech', 'all')
    args.analyse_signal_kwargs = struct();
end

args.analyse_signal_kwargs = namedargs2cell(args.analyse_signal_kwargs);

if isempty(runopt.path.oc_mna_analysis)
    analysis_folder = fullfile(runopt.path.oc_mna,'analysis');
else
    analysis_folder = runopt.path.oc_mna_analysis;
end

if ~exist(analysis_folder, 'dir')
    mkdir(analysis_folder);
end

if isempty(args.variables_to_analyze.voltage)
    voltage_to_process = {};
    
elseif ischar(args.variables_to_analyze.voltage) ...
    && strcmp(args.variables_to_analyze.voltage, 'all')

    % tmp = whos(resFiles.fVolt);
    % voltage_to_process = [ ...
    %     {tmp.name}, ...
    %     {'IHC', 'OHC'}, ...
    %     {'IHC_ic', 'IHC_ec', 'OHC_ic', 'OHC_ec', 'ScalaMedia', 'SpiralLimbus'}];
    
    voltage_to_process = [{'IHC', 'OHC'}];
else
    voltage_to_process = args.variables_to_analyze.voltage;
end

voltages = prepareVoltageStructure(voltage_to_process, 'absolute', runopt.analysisStart.ms);
voltages_relative = prepareVoltageStructure(voltage_to_process, 'relative', runopt.analysisStart.ms);

if isempty(args.variables_to_analyze.current)
    current_to_process = {};
    
elseif ischar(args.variables_to_analyze.current) ...
    && strcmp(args.variables_to_analyze.current, 'all')

    % tmp = whos(resFiles.fCurr);
    % current_to_process = {tmp.name};
    
    current_to_process = [{'IHC', 'OHC'}];
else
    current_to_process = args.variables_to_analyze.current;
end

currents = prepareCurrentStructure(current_to_process, 'absolute', runopt.analysisStart.ms);
currents_relative = prepareCurrentStructure(current_to_process, 'relative', runopt.analysisStart.ms);

%%

possible_mech_variables = struct( ...
    'BMx', struct( ...
        'id', 'BMx', ...
        'type', 'displacement', ...
        'name', 'BM displacement', ...
        'math', '\xi', ...
        'unit', 'nm'), ...
    'BMv', struct( ...
        'id', 'BMv', ...
        'type', 'velocity', ...
        'name', 'BM velocity', ...
        'math', 'd\xi/dt', ...
        'unit', 'nm/s'), ...
    'TMx', struct( ...
        'id', 'TMx', ...
        'type', 'displacement', ...
        'name', 'Stereocilia displacement', ...
        'math', '\eta', ...
        'unit', 'nm'), ...
    'TMv', struct( ...
        'id', 'TMv', ...
        'type', 'velocity', ...
        'name', 'Stereocilia velocity', ...
        'math', 'd\eta/dt', ...
        'unit', 'nm/s'));

if isempty(args.variables_to_analyze.mech)
    mech_vars_selected = {};
    
else
    if ischar(args.variables_to_analyze.mech) ...
        && strcmp(args.variables_to_analyze.mech, 'all')

        args.variables_to_analyze.mech = fieldnames(possible_mech_variables);
    end
    
    mech_vars_selected = cellfun(@(x) possible_mech_variables.(x), args.variables_to_analyze.mech);
end

if numel(mech_vars_selected) == 0
    mech_vars = prepareEmptyVariableStructure();
else
    for i = 1:numel(mech_vars_selected)
        tmp = mech_vars_selected(i);
        mech_vars(i) = prepareVariableStructure( ...
            tmp.id, 'absolute', tmp.type, tmp.unit, tmp.math, tmp.unit, 1, runopt.analysisStart.ms);
    end
end
% mech_vars(1) = prepareVariableStructure( ...
%     'BMx', 'absolute', 'mech', 'nm', '\xi', 'nm', 1, runopt.analysisStart.ms);
% mech_vars(end+1) = prepareVariableStructure( ...
%     'BMv', 'absolute', 'mech', 'nm', 'd\xi/dt', 'nm/s', 1, runopt.analysisStart.ms);
% mech_vars(end+1) = prepareVariableStructure( ...
%     'TMx', 'absolute', 'mech', 'nm', '\eta', 'nm', 1, runopt.analysisStart.ms);
% mech_vars(end+1) = prepareVariableStructure( ...
%     'TMv', 'absolute', 'mech', 'nm', 'd\eta/dt', 'nm/s', 1, runopt.analysisStart.ms);

%%

variables = [voltages, currents, voltages_relative, currents_relative, mech_vars];
all_analysis_files_exist = all([variables.analysis_file_found]);

runopt.analysis_results.oc_mna = variables;

%% Analyze results

% runopt.do.oc_oc_mna_analysis = false;

if args.all_analysis_files_must_exist
    assert(all_analysis_files_exist)
end

if  runopt.do.oc_mna_analysis && (runopt.recalculate.oc_mna || ...
        ( ~all_analysis_files_exist || runopt.recalculate.oc_mna_statistics ))

    if runopt.no_create_only_load && ~runopt.found.oc_mna
        err_str = 'Result files were not found but analysis is requested';
        if runopt.no_create_only_load == true
            warning(err_str);
            OC_results = [];
            return
        end
        error(err_str);
    end
    
    if isa(stimulus, 'PureTone')
        t1 = stimulus.main_tspan(1);
        t2 = stimulus.main_tspan(2);
    else
        t1 = topt.total.t0;
        t2 = topt.total.tf - Time(1e-10, 'ms');
        % t1 = find(runopt.analysisStart <= Time(resFiles.fT.T, resFiles.fT.unit), 1);
        % t2 = length(resFiles.fT.T);
    end
    
    t1 = find(t1 <= Time(resFiles.fT.T, resFiles.fT.unit), 1);
    t2 = find(t2 <= Time(resFiles.fT.T, resFiles.fT.unit), 1);

    if isempty(t1)
        error('runopt.analysisStart time is to long for current simulation')
    end

    clear stats
    stats.time = resFiles.fT.T(t1:t2,:);
    
    for i = 1:numel(voltage_to_process)
        
        var = voltages(i).id;
        
        sim_to_mV = Unit.conversionConstant(mnaopt.simulation_units.voltage, 'mV');
        
        TMP = sim_to_mV * getVoltage(resFiles.fVolt, var, mnaopt, 'transformed');
        VAR = TMP(t1:t2,:);
        REF = TMP(1,:); % refference (steady state)
        
        switch var
            case cellfun(@(x) sprintf('long%s', x), {'RSV','RSL','ROC'}, 'UniformOutput', false)
                xx = mnaopt.xgrid(1:end-1) + 0.5*diff(mnaopt.xgrid); % middle point
            otherwise
                xx = mnaopt.xgrid;
        end
        
        
        analyzeSignal( ...
            VAR, REF, xx, voltages(i).analysis_file, stimulus, mnaopt.samplingFrequency, mechopt, ...
            args.analyse_signal_kwargs{:});
        
        REL = VAR - REF; % relative value
        analyzeSignal( ...
            REL, zeros(size(REF)), xx, voltages_relative(i).analysis_file, stimulus, mnaopt.samplingFrequency, mechopt, ...
            args.analyse_signal_kwargs{:});
        
        % store time steps used in analysis
        save(voltages(i).analysis_file, '-struct', 'stats', '-append', '-nocompression')
        save(voltages_relative(i).analysis_file, '-struct', 'stats', '-append', '-nocompression')
        
    end
    
    for i = 1:numel(current_to_process)
        
        var = currents(i).id;
        
        VAR = resFiles.fCurr.(var);
        VAR = VAR(t1:t2,:);
        
        REF = resFiles.fCurr.(var)(1,:); % refference (steady state)
        
        analyzeSignal( ...
            VAR, REF, mnaopt.xgrid, currents(i).analysis_file, stimulus, mnaopt.samplingFrequency, mechopt, ...
            args.analyse_signal_kwargs{:});
        
        REL = VAR - REF; % relative value
        analyzeSignal( ...
            REL, zeros(size(REF)), mnaopt.xgrid, currents_relative(i).analysis_file, stimulus, mnaopt.samplingFrequency, mechopt, ...
            args.analyse_signal_kwargs{:});
        
        % store time steps used in analysis
        save(currents(i).analysis_file, '-struct', 'stats', '-append', '-nocompression')
        save(currents_relative(i).analysis_file, '-struct', 'stats', '-append', '-nocompression')
    end
    
    for i = 1:numel(mech_vars)
        
        var = mech_vars(i).id;
        
        VAR = resFiles.fMech.(var);
        VAR = VAR(t1:t2,:);
        
        % currently, we don't store relative values of MECH_VARS
        % REF = resFiles.fMech.(var)(1,:); % refference (steady state)
        
        analyzeSignal( ...
            VAR, zeros(1, size(VAR,2)), mnaopt.xgrid, mech_vars(i).analysis_file, stimulus, mnaopt.samplingFrequency, mechopt, ...
            args.analyse_signal_kwargs{:});
        
        % REL = VAR - REF; % relative value
        % analyzeSignal( ...
        %     REL, zeros(size(REF)), mnaopt.xgrid, mech_vars_relative(i).analysis_file, stimulus, mnaopt.samplingFrequency, mechopt, ...
        %     args.analyse_signal_kwargs{:});
        
        % store time steps used in analysis
        save(mech_vars(i).analysis_file, '-struct', 'stats', '-append', '-nocompression')
        % save(mech_vars_relative(i).analysis_file, '-struct', 'stats', '-append', '-nocompression')
    end
end
%%

    function voltages = prepareVoltageStructure(vars, spec, analysis_start)
        if numel(vars) == 0
            voltages = prepareEmptyVariableStructure();
        end
        for j = numel(vars):-1:1
            voltages(j) = prepareVariableStructure(vars{j}, spec, 'voltage', 'mV', 'U', 'mV', 1, analysis_start);
        end
    end
    function currents = prepareCurrentStructure(vars, spec, analysis_start)
        if numel(vars) == 0
            currents = prepareEmptyVariableStructure();
        end
        for j = numel(vars):-1:1
            currents(j) = prepareVariableStructure(vars{j}, spec, 'current', 'A', 'I', 'uA', 1e6, analysis_start);
        end
    end
    function variable = prepareVariableStructure(var, spec, type, unit, math, plot_unit, plot_unit_conversion_factor, analysis_start)
            
            if isempty(spec)
                fullid = var;
            else
                fullid = [var, spec];
            end
            
            defining_properties = struct( ...
                'id', var, ...
                'type', type, ...
                'specifier', spec, ...
                'unit', unit, ...
                'analysis_start', analysis_start);
            
            hash = DataHash(defining_properties, struct('Method', 'SHA-1', 'Format','base32'));
            
            fileid = sprintf('%s-%s-%s', type, fullid, hash);
            
            analysis_file = fullfile(analysis_folder, sprintf('%s.mat', fileid));
            
            variable = struct( ...
                'id', var, ...
                'fullid', fullid, ...
                'fileid', fileid, ...
                'type', type, ...
                'specifier', spec, ...
                'name', sprintf('%s %s %s', spec, var, type), ...
                'unit', unit, ...
                'math', math, ...
                'plot_unit', plot_unit, ...
                'plot_unit_conversion_factor', plot_unit_conversion_factor, ...
                'analysis_file', analysis_file, ...
                'analysis_file_found', exist(analysis_file, 'file'));
    end
    function variable = prepareEmptyVariableStructure()
            variable = struct( ...
                'id', [], ...
                'fullid', [], ...
                'fileid', [], ...
                'type', [], ...
                'specifier', [], ...
                'name', [], ...
                'unit', [], ...
                'math', [], ...
                'plot_unit', [], ...
                'plot_unit_conversion_factor', [], ...
                'analysis_file', [], ...
                'analysis_file_found', []);
            variable = variable([]);
    end
end
