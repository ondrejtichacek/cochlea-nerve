function [ Resistor, Capacitor, Vsource, Apical_R, Isource, numRes, numCap, numV, numI, numNode ] = circuit_elements( f_corti )
%CIRCUIT_ELEMENTS Generates arrays of circuit elements from the circuit input list represented by .txt file
%
%

%% Initialize

numRes = 0;   % Number of resistors.
numCap = 0;   % Number of capacitors.
numV = 0;     % Number of independent voltage sources
numI = 0;     % Number of independent current sources
numNode = 0;  % Number of nodes, not including ground (node 0).

Resistor =  struct('Name', {}, 'Node1', {}, 'Node2', {}, 'Value', {}, 'ValueBase', {}, 'ValueApex', {}, 'longarg', {}, 'LongValue', {});
Apical_R =  struct('Name', {}, 'Node1', {}, 'Node2', {}, 'Value', {}, 'ValueBase', {}, 'ValueApex', {}, 'longarg', {}, 'LongValue', {});
Capacitor = struct('Name', {}, 'Node1', {}, 'Node2', {}, 'Value', {}, 'ValueBase', {}, 'ValueApex', {});
Vsource =   struct('Name', {}, 'Node1', {}, 'Node2', {}, 'Value', {});
Isource =   struct('Name', {}, 'Node1', {}, 'Node2', {}, 'Value', {});

%% Read from file

% [Name, N1, N2, arg3, arg4, arg5, arg6] = textread(f_corti,'%s %s %s %s %s %s %s');
fid = fopen(f_corti);
C = textscan(fid, '%s %f %f %f %f %f %f %f %s', 'EmptyValue', NaN, 'CommentStyle', '//', 'Delimiter', {',',';'}, 'MultipleDelimsAsOne', false);
fclose(fid);

[Name, N1, N2, Value, LongArg, LongValue, ValueBase, ValueApex, SF] = C{:};

ScaleFactor = ones(size(Value));
for j = 1:length(SF)
    if ~isempty(SF{j})
        ScaleFactor(j) = eval(SF{j});
    end
end
    

%%

for i = 1:length(Name)
    switch (Name{i}(1))
        case 'R'
            numRes = numRes + 1;
            Resistor(numRes) = parseCircuitElement(Name{i}, N1(i), N2(i), Value(i), ValueBase(i), ValueApex(i), ScaleFactor(i), LongArg(i), LongValue(i));
            
            if strcmp(Name{i},'RAI')
                Apical_R(1) = Resistor(numRes);
            elseif strcmp(Name{i},'RAO')
                Apical_R(2) = Resistor(numRes);
            end
            
        case 'C'
            numCap = numCap + 1;
            Capacitor(numCap) = parseCircuitElement(Name{i}, N1(i), N2(i), Value(i), ValueBase(i), ValueApex(i), ScaleFactor(i));
            
        case 'V'
            numV = numV + 1;
            Vsource(numV) = parseCircuitElement(Name{i}, N1(i), N2(i), Value(i));
            
        case 'I'
            numI = numI + 1;
            Isource(numI) = parseCircuitElement(Name{i}, N1(i), N2(i), Value(i));
    end
    
end

numNode = max([N1;N2]);

%% parseCircuitElement function

    function Element = parseCircuitElement(Name, N1, N2, Value, ValueBase, ValueApex, ScaleFactor, LongArg, LongValue )
        
        
        if nargin < 7
            ScaleFactor = 1;
        end
        
        Element.Name = Name;
        Element.Node1 = N1;
        Element.Node2 = N2;
        Element.Value = ScaleFactor * Value;
        
        if nargin >= 6
            Element.ValueBase = ScaleFactor * ValueBase;
            Element.ValueApex = ScaleFactor * ValueApex;
        end
        if nargin >= 9
            Element.longarg = LongArg; % connective points in the sparse matrix must be atributed to resistor components
            Element.LongValue = LongValue;
        end
    end

end