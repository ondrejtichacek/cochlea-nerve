function [ Voltage, sel, f, operand, fun ] = getVoltage( y, NAME, mnaopt, flag, tpoints )
%GETVOLTAGE 
arguments
    y
    NAME
    mnaopt
    flag = 'original'
    tpoints = ':'
end

numNode = mnaopt.circuit.num_node;
numV = mnaopt.circuit.num_v;
Numstacks = mnaopt.Numstacks;

operand = '';
fun = [];

if isempty(NAME)
    Voltage = zeros(size(y,1), Numstacks, numNode);
    sel{numNode} = [];
    for jj = 1:numNode
        sel{jj} = ((1:Numstacks)-1)*(numNode+numV)+jj;
        Voltage(:, 1:Numstacks, jj) = y(:,sel{jj});
    end
else
    
    ele = mnaopt.circuit.select('type resistor and special longitudinal');
    long_resistors = {ele.name};
    
    switch NAME
        case cellfun(@(x) sprintf('N%d', x), num2cell(1:numNode), 'UniformOutput', false)
            jj = str2double(NAME(2:end));
            [Voltage, sel, f] = getNodeVoltage(y, jj);
       
        case long_resistors
            jj = mnaopt.circuit.get_element_by_name(NAME).node1;
            
            [tmp, sel, f] = getNodeVoltage(y, jj);
            
            sel1 = sel(1:end-1);
            sel2 = sel(2:end);
            
            Voltage(:,1:Numstacks-1) = diff(tmp,1,2);
            sel = {sel2, sel1};
            operand = '-';
            
        case {mnaopt.circuit.resistors.name}
            
            if mnaopt.circuit.get_element_by_name(NAME).special == "longitudinal"
            else
                j1 = mnaopt.circuit.get_element_by_name(NAME).node1;
                j2 = mnaopt.circuit.get_element_by_name(NAME).node2;

                [Voltage, sel, f, operand, fun] = nodeDiff(y, j1, j2);
            end
            
        case 'ScalaMedia'
            jj = 1;
            [Voltage, sel, f] = getNodeVoltage(y, jj);
            
        case 'IHC_ic'
            jj = mnaopt.circuit.get_element_by_name('CAI').node2;
            [Voltage, sel, f] = getNodeVoltage(y, jj);
            
        case 'OHC_ic'
            jj = mnaopt.circuit.get_element_by_name('CAO').node2;
            [Voltage, sel, f] = getNodeVoltage(y, jj);
            
        case 'IHC_ec'
            jj = mnaopt.circuit.get_element_by_name('CBI').node2;
            [Voltage, sel, f] = getNodeVoltage(y, jj);
            
        case 'OHC_ec'
            jj = mnaopt.circuit.get_element_by_name('CBO').node2;
            [Voltage, sel, f] = getNodeVoltage(y, jj);
            
        case 'SpiralLigament'
            jj = 8;
            [Voltage, sel, f] = getNodeVoltage(y, jj);
            
        case 'IHC'
            %j1 = mnaopt.circuit.get_element_by_name('RAI').node1;
            %j2 = mnaopt.circuit.get_element_by_name('RBI').node2;
            
            %j1 = mnaopt.circuit.get_element_by_name('RAI').node2;
            %j2 = mnaopt.circuit.get_element_by_name('VI').node2;
            
            j1 = mnaopt.circuit.get_element_by_name('CBI').node1;
            j2 = mnaopt.circuit.get_element_by_name('CBI').node2;
            
            [Voltage, sel, f, operand, fun] = nodeDiff(y, j1, j2);
            % [Voltage, sel, f, operand, fun] = nodePlus(y, j1, j2);
            
%             Voltage = Voltage - 0.07;
%             fun = @(y) fun(y) - 0.07;
            
        case 'OHC'
            %j1 = mnaopt.circuit.get_element_by_name('RAO').node1;
            %j2 = mnaopt.circuit.get_element_by_name('RBO').node2;
            
            %j1 = mnaopt.circuit.get_element_by_name('RAO').node2;
            %j2 = mnaopt.circuit.get_element_by_name('VO').node2;
            
            j1 = mnaopt.circuit.get_element_by_name('CBO').node1;
            j2 = mnaopt.circuit.get_element_by_name('CBO').node2;
            
            [Voltage, sel, f, operand, fun] = nodeDiff(y, j1, j2);
            % [Voltage, sel, f, operand, fun] = nodePlus(y, j1, j2);
            
%             Voltage = Voltage - 0.07;
%             fun = @(y) fun(y) - 0.07;

        case 'IHC_MET'
            
            j1 = mnaopt.circuit.get_element_by_name('CAI').node1;
            j2 = mnaopt.circuit.get_element_by_name('CAI').node2;
            
            [Voltage, sel, f, operand, fun] = nodeDiff(y, j1, j2);
            
        case 'OHC_MET'
            
            j1 = mnaopt.circuit.get_element_by_name('CAO').node1;
            j2 = mnaopt.circuit.get_element_by_name('CAO').node2;
            
            [Voltage, sel, f, operand, fun] = nodeDiff(y, j1, j2);

        otherwise
            error('Unknown option %s', NAME);
    end
end
    function [Voltage, sel, f] = getNodeVoltage(y, index)
        switch flag
            case 'original'
                sel = ((1:Numstacks)-1)*(numNode+numV) + index;
                if isempty(y)
                    Voltage = [];
                else
                    Voltage(:,1:Numstacks) = y(:,sel);
                end
                f = '';
            case 'transformed'
                sel = 1:Numstacks;
                f = sprintf('N%d', index);
                if isempty(y)
                    Voltage = [];
                else
                    Voltage = y.(f)(tpoints,:);
                end
            otherwise
                error('Unknown flag %s');
        end
    end
    function [Voltage, sel, f, operand, fun] = nodeDiff(y, j1, j2)
        operand = '-';
        
        if j1 == 0
            [V2, s2] = getNodeVoltage(y, j2);
            
            Voltage = -V2;
            sel = s2;
            
            fun = @(V2) -V2;
            
            f = sprintf('N%d', j2);
        elseif j2 == 0
            [Voltage, sel] = getNodeVoltage(y, j1);
            
            fun = @(V1) V1;
            
            f = sprintf('N%d', j1);
        else
            [V1, s1] = getNodeVoltage(y, j1);
            [V2, s2] = getNodeVoltage(y, j2);

            Voltage = V1 - V2;
            sel = {s1, s2};
            
            fun = @(V1,V2) V1 - V2;
            
            f = {sprintf('N%d', j1), sprintf('N%d', j2)};
        end
    end
    function [Voltage, sel, f, operand, fun] = nodePlus(y, j1, j2)
        operand = '+';
        
        if j1 == 0
            [Voltage, sel] = getNodeVoltage(y, j2);
            
            fun = @(V2) V2;
            
            f = sprintf('N%d', j2);
        elseif j2 == 0
            [Voltage, sel] = getNodeVoltage(y, j1);
            
            fun = @(V1) V1;
            
            f = sprintf('N%d', j1);
        else
            [V1, s1] = getNodeVoltage(y, j1);
            [V2, s2] = getNodeVoltage(y, j2);
            
            Voltage = V1 + V2;
            sel = {s1, s2};
            
            fun = @(V1,V2) V1 + V2;
            
            f = {sprintf('N%d', j1), sprintf('N%d', j2)};
        end
    end

end
