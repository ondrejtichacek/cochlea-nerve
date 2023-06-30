classdef CompoundSignal < Opt %< Signal & SignalWithAmplitude
    %COMPOUNDSIGNAL
    
    properties
        compounds
        nc
    end
    properties (Hidden)
        definingProperties = { ...
            'compounds',
            };

    end
    properties (Dependent)
        requiredParameters
        
        tspan
        audio
        audiotime
        name
        amplitude
        frequency
        samplingFrequency

        t0
        tf
        fs
    end
    methods (Static)
        function par = getrequiredParameters()
            par = {'compunds'};
        end
    end
    methods
        function obj = CompoundSignal(compounds)
            
            obj.compounds = compounds;
            obj.nc = numel(compounds);
        end
        function t0 = get.t0(obj)
            t0 = obj.tspan(1);
        end
        function tf = get.tf(obj)
            tf = obj.tspan(2);
        end
        function fs = get.fs(obj)
            fs = obj.samplingFrequency;
        end
        % =================================================================
        % INTERFACES
        function par = get.requiredParameters(obj)
            par = obj.getrequiredParameters();
        end
        function intensity = relative_intensity(obj)
            intensity = 10.^([obj.compounds.spl]/20);
            intensity = intensity / max(intensity);
        end
        function val = get.tspan(obj)
            val = {obj.compounds.tspan};
            for i = 2:obj.nc
                assert(all(size(val{1}) == size(val{i})))
                assert(all(val{1} == val{i}))
            end
            val = val{1};
        end
        function val = get.audio(obj)
            val = {obj.compounds.audio};
            for i = 2:obj.nc
                assert(all(size(val{1}) == size(val{i})))
                assert(size(val{i},1) == 1)
            end

            intensity = obj.relative_intensity();

            for i = 1:obj.nc
                val{i} = val{i} * intensity(i);
            end

            val = sum(cat(1,val{:}));
        end
        function val = get.audiotime(obj)
            val = {obj.compounds.audiotime};
            for i = 2:obj.nc
                assert(all(size(val{1}) == size(val{i})))
                assert(all(val{1} == val{i}))
            end
            val = val{1};
        end
        function val = get.name(obj)
            val = [obj.compounds.name];
        end
        function val = get.amplitude(obj)
            val = max([obj.compounds.amplitude]);
        end
        function val = get.frequency(obj)
            val = [obj.compounds.frequency];
        end
        function val = get.samplingFrequency(obj)
            val = [obj.compounds.samplingFrequency];
            assert(all(val(1) == val))
            val = val(1);
        end
        function t = zeroDuration(obj)
            t = Time(+inf, 's');
            for i = 1:numel(obj.compounds)
                t = min([t, obj.compounds(i).zeroDuration]);
            end
        end
        function t = default_analysis_start_time(obj)
            t = Time(+inf, 's');
            for i = 1:numel(obj.compounds)
                t = min([t, obj.compounds(i).default_analysis_start_time]);
            end
        end
        function signal = eval(obj, time)
    
            A_rel = obj.relative_intensity();

            signal = A_rel(1) * obj.compounds(1).eval(time);
            for i = 2:numel(obj.compounds)
                signal = signal + A_rel(i) * obj.compounds(i).eval(time);
            end
        end
        function fun = eval_fcn(obj, time_unit)
            error('Need to check if amplitude scaling is done properly')
    
            A_rel = obj.relative_intensity();

            for i = 1:obj.nc
                signal_fun_partial{i} = obj.compounds(i).eval_fcn(time_unit);
            end 
            fun = @(t) cellfun(@(f) f(t), signal_fun_partial, ...
                'UniformOutput', true);
        end
        function fun = envelope_fcn(obj, time_unit)

            A_rel = obj.relative_intensity();

            % error('Need to check if amplitude scaling is done properly')
            envelopes = cell(1,obj.nc);
            for i = 1:obj.nc 
                envelopes{i} = obj.compounds(i).envelope_fcn(time_unit);
            end
            fun = @(t) obj.compound_envelope(t, envelopes, A_rel);
        end
        function plot(obj)
            
            % audiotime = obj.audiotime;

            sig = VarSignal( ...
                'originalAudio', obj.audio, ...
                'originalSamplingFrequency', obj.samplingFrequency, ...
                'spl', max([obj.compounds.spl]), ...
                'samplingFrequency', obj.samplingFrequency, ...
                'tspan', obj.tspan, ...
                'name', 'tmp');

            plot(sig)
        end

        function play(obj)
            
            % audiotime = obj.audiotime;

            sig = VarSignal( ...
                'originalAudio', obj.audio, ...
                'originalSamplingFrequency', obj.samplingFrequency, ...
                'spl', max([obj.compounds.spl]), ...
                'samplingFrequency', obj.samplingFrequency, ...
                'tspan', obj.tspan, ...
                'name', 'tmp');

            play(sig)
        end
    end

    methods (Static)
        function e = compound_envelope(t, envelopes, A_rel)
            e = zeros(size(t));
            for i = 1:numel(envelopes)
                e = max(e, envelopes{i}(t) * A_rel(i));
            end
        end
    end
    
end
