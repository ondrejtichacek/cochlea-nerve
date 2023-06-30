classdef Circuit < Opt
    %CIRCUIT 
    
    properties
        elements
    end
    properties (Dependent)
        num_res
        num_ind
        num_cap
        num_v
        num_i
        num_node
        resistors
        inductors
        capacitors
        vsources
        isources
    end
    properties (Hidden)
        definingProperties = { ...
            'elements', ...
            }
    end
    properties (Constant, Hidden)
        requiredParameters = {};
    end
    methods        
        function val = get_element_by_name(obj, name)
            val = obj.elements(find_elements(obj, 'name', name));
        end
        function add_element(obj, element)
            if isempty(obj.elements)
                obj.elements = element;
            else
                obj.elements(end+1) = element;
            end
        end
        function count = count_elements(obj, prop, val)
            if ischar(val)
                count = sum(strcmp(val, {obj.elements.(prop)}));
            else
                count = sum(val == [obj.elements.(prop)]);
            end
        end
        function ind = find_elements(obj, prop, val)
            if ischar(val)
                ind = strcmp(val, {obj.elements.(prop)});
            else
                ind = val == [obj.elements.(prop)];
            end
        end
        function val = select(obj, selection_string)
            parts = strsplit(selection_string, ' and ');
            
            f = @(A,B) and(A,B);

            ind = ones(size(obj.elements));
            
            for i = 1:numel(parts)
                p = strsplit(parts{i});
                
                if strcmp(p{2}, 'not')
                    p = p([1,3]);
                    ff = @(A,B) f(A,not(B));
                else
                    ff = f;
                end
                
                assert(numel(p) == 2)
                
                ii = obj.find_elements(p{:});
                
                ind = ff(ind, ii);
            end
            val = obj.elements(ind);
        end
        
        function val = get.num_node(obj)
            val = max([obj.elements.node1, obj.elements.node2]);
        end
        function val = get.num_res(obj)
            val = obj.count_elements('type', 'resistor');
        end
        function val = get.num_ind(obj)
            val = obj.count_elements('type', 'inductor');
        end
        function val = get.num_cap(obj)
            val = obj.count_elements('type', 'capacitor');
        end
        function val = get.num_v(obj)
            val = obj.count_elements('type', 'vsource');
        end
        function val = get.num_i(obj)
            val = obj.count_elements('type', 'isource');
        end
        function val = get.resistors(obj)
            val = obj.elements(obj.find_elements('type', 'resistor'));
        end
        function val = get.inductors(obj)
            val = obj.elements(obj.find_elements('type', 'inductor'));
        end
        function val = get.capacitors(obj)
            val = obj.elements(obj.find_elements('type', 'capacitor'));
        end
        function val = get.vsources(obj)
            val = obj.elements(obj.find_elements('type', 'vsource'));
        end
        function val = get.isources(obj)
            val = obj.elements(obj.find_elements('type', 'isource'));
        end
    end
    
end

