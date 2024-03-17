classdef SignalWithAmplitude < Opt

    properties
        amplitude
    end
    properties (Dependent)
        frequency
        spl
        phon
    end
    properties (Access=private)
        frequency_
        spl_
        phon_
    end

    methods
        % =================================================================
        function obj = SignalWithAmplitude()
        end
        % =================================================================
        % GET FUNCTIONS
        % -----------------------------------------------------------------
        function val = get.frequency(obj)
            val = obj.frequency_;
        end
        function val = get.phon(obj)
            val = obj.phon_;
        end
        function val = get.spl(obj)
            val = obj.spl_;
        end
        % =================================================================
        % SET FUNCTIONS
        % -----------------------------------------------------------------
        function set.frequency(obj, val)
            if ~isempty(val)
                Opt.myassert(val, 'scalar', 'Frequency')
            
                obj.frequency_ = val;
                if ~isempty(obj.phon)
                    obj.amplitude = iosr.auditory.iso226(obj.phon, val.Hz);
                end
            end
        end
        function set.phon(obj, val)
            if ~isempty(val)
                Opt.myassert(val, 'scalar', 'numeric')

                obj.phon_ = val;
                if ~isempty(obj.frequency)
                    obj.amplitude = iosr.auditory.iso226(val, obj.frequency.Hz);
                end
            end
        end
        function set.spl(obj, val)
            if ~isempty(val)
                Opt.myassert(val, 'scalar', 'numeric')

                obj.spl_ = val;
                obj.amplitude = val;
            end
        end
    end
end
