classdef AnalysisFileWrapper < matlab.mixin.Copyable
    %ANALYSISFILEWRAPPER 
    
    properties       
        analysis_files
        shape
    end
    
    methods
        function obj = AnalysisFileWrapper(shape)
            %ANALYSISFILEWRAPPER 
            obj.shape = shape;
            obj.analysis_files = cell(obj.shape{:});
        end
        
        function add_analysis_file(obj, mf, ind, args)
            arguments
                obj
                mf
                ind
                args.load = false
            end
            if args.load == true
                mf = load(mf.Properties.Source);
            end
            obj.analysis_files{ind{:}} = mf;            
        end
        function list_props(obj)
            whos(obj.analysis_files{1})
        end
        function val = get_prop(obj, VAR, array_type, varargin)
            if isempty(varargin)
                switch array_type
                    case 'array'
                        val = zeros(obj.shape{:});
                        for i = 1:numel(obj.analysis_files)                            
                            if isa(obj.analysis_files{i}, 'matlab.io.MatFile') ...
                                    || isa(obj.analysis_files{i}, 'struct')
                                val(i) = obj.analysis_files{i}.(VAR);
                            elseif isnan(obj.analysis_files{i})
                                val(i) = NaN;
                            else
                                error('unsupported class');
                            end
                        end
                    case 'cell'
                        val = cell(obj.shape{:});
                        for i = 1:numel(obj.analysis_files)
                            if isa(obj.analysis_files{i}, 'matlab.io.MatFile') ...
                                    || isa(obj.analysis_files{i}, 'struct')
                                val{i} = obj.analysis_files{i}.(VAR);
                            elseif isnan(obj.analysis_files{i})
                                val{i} = NaN;
                            else
                                error('unsupported class');
                            end
                        end
                    otherwise
                        error('Unrecognised type %s', array_type)
                end
            else
                V = obj.analysis_files{varargin{:}};
                if isa(V, 'matlab.io.MatFile') || isstruct(V)
                    val = V.(VAR);
                elseif isnan(V)
                    val = NaN;
                else
                    error('unsupported class');
                end
            end
        end
    end
end

