function [ stimulus_out ] = MiddleEar( ...
    stimulus_in, middle_ear_filter, middle_ear_params, ...
    outer_ear_params, ...
    tspanopt, topt, midopt, runopt, opt, memopt, paropt )
%MIDDLEEAR
arguments
    stimulus_in
    middle_ear_filter (1,:) char
    middle_ear_params (1,1) struct
    outer_ear_params (1,1) struct
    tspanopt
    topt (1,1) timeOpt
    midopt (1,1) midOpt
    runopt (1,1) runOpt
    opt (1,1) globalOpt
    memopt (1,1) memOpt
    paropt (1,1) parOpt
end

num_signal = numel(stimulus_in);

fs = [stimulus_in.samplingFrequency];

assert(all(fs == fs(1)), 'Sampling frequency must be the same');

fs = fs(1);

switch middle_ear_filter
    case 'none'
        
        A = stimulus_in.amplitude;

        % we need two extra samples over the time span to compensate for 
        % the second numerical derivative.
        extra_points = 2;
        
        audiotime = tspanopt.generateTimeSamples(fs, extra_points);

        audio = zeros(size(audiotime));

        for i = 1:num_signal
            audio = audio + stimulus_in(i).eval(audiotime);
        end
        
        AM = mech.db2input(A, ...
            middle_ear_params.AMo,...
            middle_ear_params.L); % [length_unit * s^2]
        
        if isa(stimulus_in, 'Click2')
            acceleration = - AM * audio * (stimulus_in.om.Hz^2);
            accelerationtime = audiotime;
        else
            acceleration = AM * diff(audio,2) * fs.Hz^2;
            % accelerationtime = audiotime(1:end-2);
            acceleration = [0, 0, acceleration];
            accelerationtime = audiotime;
        end
        stimulus_out = StapesSignal( ...
            'audio', audio, ...
            'audiotime', audiotime, ...
            'accelerationtime', accelerationtime, ...
            'samplingFrequency', fs, ...
            'acceleration', acceleration);
        
    case {'lopezpoveda', 'jepsen'}
        
        if num_signal > 1
            error('Not yet implemented')
        end

        A = stimulus_in.amplitude;

        % we need two extra samples over the time span to compensate for 
        % the first numerical derivative.
        extra_points = 1;
        
        audiotime = tspanopt.generateTimeSamples(fs, extra_points);
        audio = stimulus_in.eval(audiotime);
        
        ME = MiddleEarFilter( ...
            'samplingFrequency', fs, ...
            'middle_ear_filter', middle_ear_filter);
        
        [audio, velocity, acceleration] = ME.applyFilter(audio, A, fs);
        
        % accelerationtime = audiotime(1:end-1);
        accelerationtime = audiotime(1:end-2);
        
        
        stimulus_out = StapesSignal( ...
            'audio', audio, ...
            'velocity', velocity, ...
            'accelerationtime', accelerationtime, ...
            'samplingFrequency', fs, ...
            'acceleration', acceleration);

    case {'PBLL', 'Zwislocki'}
        
%         if num_signal > 1
%             error('Not yet implemented')
%         end

        resFiles = ResCacheMiddleEar( ...
            stimulus_in, topt, midopt, runopt, opt, memopt, paropt);
        
        % The result of ME that we care about is the voltage equivalent to
        % the displacement of the stapes (oval window)
        V_OW = resFiles.fy.y(:,midopt.IND) - resFiles.fy.y(1,midopt.IND);
        t = resFiles.ft.t;
        
        F = griddedInterpolant(t, V_OW);
        
        % we need two extra samples over the time span to compensate for 
        % the second numerical derivative.
        extra_points = 2;
        
        stapes_time = tspanopt.generateTimeSamples(fs, extra_points);
        stapes_signal = F(stapes_time.(midopt.simulation_units.time));        
        
        k = Unit.conversionConstant(midopt.simulation_units.voltage, 'V');
        
        stapes_signal = stapes_signal * (k * middle_ear_params.vd_OW); % scaling factor of oval window
        
        % stapes_signal_meaning = 'displacement';
        stapes_signal_meaning = 'velocity';
        % stapes_signal_meaning = 'acceleration';
        
        
        switch stapes_signal_meaning 
            case 'displacement'
                acceleration = diff(stapes_signal,2) * fs.Hz^2;
                accelerationtime = stapes_time(1:end-2);
            case 'velocity'
                acceleration = diff(stapes_signal) * fs.Hz;
                accelerationtime = stapes_time(1:end-1);
            case 'acceleration'
                acceleration = stapes_signal;
                accelerationtime = stapes_time;
            otherwise
                error('undefined value of stapes_signal_meaning: %s', stapes_signal_meaning);
        end
        
        oe = outer_ear_params;
        
        % we can apply the filter here since both ME and OE are linear
        switch oe.identifier
            case 'Meddis'
                kwargs = namedargs2cell(oe.args);
                [acceleration] = OuterEar_Meddis(acceleration, fs.Hz, kwargs{:});
            case 'none'
            otherwise
                error('Unknown Outer Ear identifier %s.', oe.identifier)
        end
        
        stimulus_out = StapesSignal( ...
            'accelerationtime', accelerationtime, ...
            'samplingFrequency', fs, ...
            'acceleration', acceleration);
        
    otherwise 
        error('unknown middle_ear_filter %s', middle_ear_filter);
end

end
