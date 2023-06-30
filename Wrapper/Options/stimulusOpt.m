classdef stimulusOpt < Opt

    properties
        contTime
        samplingFrequency
    end
    properties (Dependent, Hidden)
        requiredParameters
    end
    methods (Static)
        function par = getrequiredParameters()
            par = {};
        end
    end
    methods
        function par = get.requiredParameters(obj)
            par = obj.getrequiredParameters();
        end
    end
    methods (Static)
        function stimtype = parseStimType(varargin)

            stimtype = '';

            for i = (length(varargin)-1) : -1 : 1
                if isa(varargin{i}, 'char') && strcmp(varargin{i}, 'type')
                    stimtype = varargin{i+1};
                    break;
                end
            end

            if isempty(stimtype)
                for i = length(varargin) : -1 : 1
                    if isa(varargin{i}, 'stimulusOpt')
                        stimtype = varargin{i}.type;
                        break;
                    end
                end
            end

            if isempty(stimtype)
                error('Stimulus type must be specified')
            end
        end
    end

    methods
        % -----------------------------------------------------------------
        function obj = stimulusOpt(varargin)
            % gather name-value pair input arguments
            opt = struct(varargin{:});

            % check if all required properties are set
            obj.checkRequired(opt);

            % assign the struct values to the object
            obj.expandProperties(opt);
        end
    end

end
