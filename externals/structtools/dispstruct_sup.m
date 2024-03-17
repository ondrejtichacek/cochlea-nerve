classdef dispstruct_sup < matlab.mixin.CustomDisplay
    properties
        S
    end
    methods
    function obj = dispstruct_sup(S)
        obj.S = S;
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

