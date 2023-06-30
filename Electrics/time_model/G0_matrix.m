function [ G ] = G0_matrix( Resistor, numRes, numNode, ...
        NumstacksScaling, dependencies )
%G0_MATRIX calculates matrix of R conductance (not time-dependent)

G = zeros(numNode,numNode);

%% Fill the G matrix with transducer resistros from netlist

allowed_deps = {'xpos', 'vihc_ss', 'vohc_ss'};

for i = 1:numRes
    
    if strcmp(Resistor(i).special, 'longitudinal')
        continue
    end
    
    % disp(Resistor(i))
    
    deps_ok = true;
    % allow only resistors that depend on nothing or xpos
    for j = 1:numel(Resistor(i).dep)
%         if all(strcmp(Resistor(i).dep(j), allowed_deps) == 0)
%             deps_ok = false;
%         end
        if is_in(Resistor(i).dep{j}, allowed_deps)
        else
            deps_ok = false;
        end
    end
    
    if deps_ok == false
        continue
    end
    
    n1 = Resistor(i).node1;
    n2 = Resistor(i).node2;

    dep = namedargs2cell(dependencies);

    R = Resistor(i).value(dep{:});

    Resistor(i).debug.value_assigned = 'G0';
    
    switch Resistor(i).si_unit
        case 'S'
            g = R;
        case 'Ohm'
            g = 1/R;
        otherwise
            error('Unknown unit %s.', Resistor(i).si_unit);
    end

    g = g * NumstacksScaling; % parallel resistors

    % If neither side of the element is connected to ground
    % then subtract it from appropriate location in matrix.
    if (n1 ~= 0) && (n2 ~= 0)
        G(n1,n2) = G(n1,n2) - g;
        G(n2,n1) = G(n2,n1) - g;
    end

    % If node 1 is not connected to ground, add element to diagonal
    % of matrix.
    if (n1 ~= 0)
        G(n1,n1) = G(n1,n1) + g;
    end
    % Ditto for node 2.
    if (n2 ~= 0)
        G(n2,n2) = G(n2,n2) + g;
    end        
    
end
