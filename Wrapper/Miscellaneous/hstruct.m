classdef hstruct < handle & dynamicprops
  properties
  end
  methods
    function obj = hstruct(name, value)
        arguments (Repeating)
            name
            value
        end
        for i = 1:numel(name)
            obj.addprop(name{i});
            obj.(name{i}) = value{i};
        end
    end
  end
end