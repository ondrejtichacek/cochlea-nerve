classdef WrapFile < Opt
    %WRAPFILE 
    
    properties
        name
        extension
        dir
        type
        handle
        save_method
    end
    properties (Dependent)
        path
    end
    properties (Hidden)
        definingProperties = {}
        requiredParameters = {}
    end
    methods
        function obj = WrapFile(varargin)
            opt = struct(varargin{:});            
            obj.checkRequired(opt);            
            obj.expandProperties(opt);
        end
        function fopen(obj, varargin)
            if ~isempty(obj.handle)
                warning('File already opened, closing handle first\n')
                obj.fclose();
            end
            obj.handle = fopen(obj.path, varargin{:});
        end
        function fclose(obj)
            fclose(obj.handle);
            obj.handle = [];
        end
        function val = fread(obj, varargin)
            val = fread(obj.handle, obj.type, varargin{:});
        end
        function fwrite(obj, val, varargin)
            fwrite(obj.handle, val, obj.type, varargin{:});
        end
        function val = get.path(obj)
            val = fullfile(obj.dir, [obj.name,'.',obj.extension]);
            if find(val == '~', 1)
                val = strrep(val, '~', getenv('HOME'));
            end
        end
    end
end

