classdef Opt < matlab.mixin.Copyable & matlab.mixin.CustomDisplay % & matlab.mixin.SetGet
    %OPT
    properties (Dependent, SetAccess = private, Hidden)
        DefiningData
    end
    properties (Abstract, Hidden)
        definingProperties
    end
    
    methods       
        function hash = uint8(obj)
            
            OPT.Method = 'SHA-1';
            OPT.Format = 'uint8';
            hash = DataHash( ...
                struct(...
                    'defining_data', [obj.DefiningData], ...
                    'metadata', struct( ...
                        'class', class(obj), ...
                        'size', size(obj))), ...
                OPT);
            
        end
        
        function obj = Opt()
        end        
                
        function update_struct(obj, p, fields, value, subfield)
            arguments
                obj
                p
                fields
                value
                subfield (1,:) char = ''
            end
           % updates the structure property:
           % p.(fields) = value
           
           if isempty(subfield)
               for i = 1:numel(fields)
                   obj.(p).(fields{i}) = value;
               end
           else
               for i = 1:numel(fields)
                   obj.(p).(subfield).(fields{i}) = value;
               end
           end
        end
        
        function update(obj,varargin)
            
%             allparams = catstruct( ...
%                 obj.requiredParameters, ...
%                 obj.optionalParameters, ...
%                 obj.defaultParameters );
%             
%             opt = optParser( ...
%                 struct(), ...
%                 allparams, ...
%                 struct(), ...
%                 varargin{:} );

            % gather name-value pair input arguments
            opt = struct(varargin{:});
            
            % assign the struct values to the object
            obj.expandProperties(opt);
            
        end
        function set(obj, names, values)
            arguments
                obj Opt
            end
            arguments (Repeating)
                names char
                values
            end
            for i = 1:numel(names)
                name = names{i};
                val = values{i};
                
                nameparts = strsplit(name, '.');
                
                if numel(nameparts) == 1
                    n1 = nameparts{1};
                    obj.(n1) = val;
                elseif numel(nameparts) == 2
                    n1 = nameparts{1};
                    n2 = nameparts{2};
                    obj.(n1).(n2) = val;
                else
                    error('We can currently do only depth-1 assignment');
                end
            end
        end
        function obj = set_return(obj, varargin)
            obj.set(varargin{:});
        end
        function updateWithObject(old, new)
            
            
            if IsOctave()
                f = fieldnames(struct(new));
            else
                f = properties(new);
            end
            
            for i = 1:numel(f)
                
                mp = findprop(old,f{i});
                if strcmp(mp.SetAccess, 'public')
                    if ~isempty(new.(f{i}))
                        old.(f{i}) = new.(f{i});
                    end
                end
            end
            
        end
        
        function expandProperties(obj,opt)
            
            f = fieldnames(opt);
            
            for i = 1:numel(f)
                obj.(f{i}) = opt.(f{i});
            end
                        
        end
        
        function checkRequired(obj, opt, par)
            if nargin < 3
                par = obj.requiredParameters;
            end
            
            % Check if any required parameter is missing
            if ~isempty(par)
                reqParamMissing = ~isfield(opt, par);

                if any(reqParamMissing)
                    reqParamMissingNames = sprintf('%s, ', par{reqParamMissing});
                    error('Required parameter `%s` not set', reqParamMissingNames(1:end-2));
                end
            end
        end
        
        function val = get.DefiningData(obj)
            
            val = struct();
            for i = 1:numel(obj.definingProperties)
                val.(obj.definingProperties{i}) = obj.(obj.definingProperties{i});
            end
            
        end
        
        function extendDefiningProperties(obj, props)
            obj.definingProperties = unique([obj.definingProperties, props]);
        end
        
    end
    methods (Static)
        function myassert(val, varargin)
            for i = 1:length(varargin)
                switch varargin{i}
                    case 'cell'
                        assert(isa(val,'cell'), ...
                            'Input must be a cell, not a `%s`.', class(val))
                    case 'char'
                        assert(isa(val,'char'), ...
                            'Input must be a char, not a `%s`.', class(val))
                    case 'logical_convertible'
                        assert(islogical(logical(val)), ...
                            'Input must be convertible to logical, not a `%s`.', class(val))
                    case 'natural'
                        assert(all(abs(val) == round(val)), ...
                            'Input must be an natural number.')
                    case 'numeric'
                        assert(isnumeric(val), ...
                            'Input must be numeric; current type `%s`.', class(val))
                    case 'positive'
                        assert(all(val > 0), ...
                            'Input must be positive.')
                    case 'scalar'
                        assert(isscalar(val), ...
                            'Input must be a scalar; current size [ %s].', sprintf('%d ', size(val)))
                    case 'tspanOpt'
                        assert(isa(val,'tspanOpt'), ...
                            'Input must be a tspanOpt, not a `%s`.', class(val))
                    case 'Time'
                        assert(isa(val,'Time'), ...
                            'Input must be a Time, not a `%s`.', class(val))
                    case 'Frequency'
                        assert(isa(val,'Frequency'), ...
                            'Input must be a Frequency, not a `%s`.', class(val))
                    otherwise
                        error('Undefined type parameter `%s`.', varargin{i})
                end
            end
        end
        
        function val = sprint_set(set_cell)
            spaces = repmat({', '}, size(set_cell));
            spaces{end-1} = ' or ';
            spaces{end} = '';
            joined = [set_cell; spaces];
            val = strcat([joined{:}]);
        end
        
    end
    methods (Access = protected)
       function propgrp = getPropertyGroups(obj)
          if ~isscalar(obj)
             propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
          else
             
             prop = properties(obj);
             s = struct();
             for i = 1:length(prop)
                 if any(strcmp(superclasses(obj.(prop{i})), 'Unit'))
                     if numel(obj.(prop{i}))  == 1
                        s.(prop{i}) = sprintf('%s [1x1 %s]', obj.(prop{i}), class(obj.(prop{i})));
                     else
                        s.(prop{i}) = obj.(prop{i});
                     end
                 else
                    s.(prop{i}) = obj.(prop{i});
                 end
             end
             
             propgrp(1) = matlab.mixin.util.PropertyGroup(s, '');
             
          end
       end
    end
    
end
