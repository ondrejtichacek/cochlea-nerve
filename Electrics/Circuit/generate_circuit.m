function [circuit] = generate_circuit( connection, ...
    resistors, inductors, capacitors, vsources, isources, ...
    units, simulation_units)
% GENERATE_CIRCUIT

circuit = Circuit();

si_units = struct( ...
    'conductance', 'S', ...
    'inductance', 'H', ...
    'capacitance', 'F', ... 
    'voltage', 'V', ...
    'current', 'A');

% -------------------------------------------------------------------------
% RESISTORS

k = Unit.conversionConstant(units.conductance, simulation_units.conductance);

for i = 1:numel(resistors)
    el = CircuitElement('type', 'resistor', resistors{i}{:});
    
    el.node1 = connection.(el.name)(1);
    el.node2 = connection.(el.name)(2);
    
    el.ref_value = el.ref_value * k;
    
    el.unit = simulation_units.conductance;
    el.si_unit = si_units.conductance;
    
    circuit.add_element(el);
end

% -------------------------------------------------------------------------
% INDUCTORS

k = Unit.conversionConstant(units.inductance, simulation_units.inductance);

for i = 1:numel(inductors)
    el = CircuitElement('type', 'inductor', inductors{i}{:});
    
    el.node1 = connection.(el.name)(1);
    el.node2 = connection.(el.name)(2);
    
    el.ref_value = el.ref_value * k;
    
    el.unit = simulation_units.inductance;
    el.si_unit = si_units.inductance;
    
    circuit.add_element(el);
end

% -------------------------------------------------------------------------
% CAPACITORS

k = Unit.conversionConstant(units.capacitance, simulation_units.capacitance);

for i = 1:numel(capacitors)
    el = CircuitElement('type', 'capacitor', capacitors{i}{:});
    
    el.node1 = connection.(el.name)(1);
    el.node2 = connection.(el.name)(2);
    
    el.ref_value = el.ref_value * k;
    
    el.unit = simulation_units.capacitance;
    el.si_unit = si_units.capacitance;
    
    circuit.add_element(el);
end

% -------------------------------------------------------------------------
% VSOURCES

k = Unit.conversionConstant(units.voltage, simulation_units.voltage);

for i = 1:numel(vsources)
    el = CircuitElement('type', 'vsource', vsources{i}{:});
    
    el.node1 = connection.(el.name)(1);
    el.node2 = connection.(el.name)(2);
    
    el.ref_value = el.ref_value * k;    
    
    el.unit = simulation_units.voltage;
    el.si_unit = si_units.voltage;
    
    circuit.add_element(el);
end

% -------------------------------------------------------------------------
% ISOURCES

k = Unit.conversionConstant(units.current, simulation_units.current);

for i = 1:numel(isources)
    el = CircuitElement('type', 'isource', isources{i}{:});
    
    el.node1 = connection.(el.name)(1);
    el.node2 = connection.(el.name)(2);
    
    el.ref_value = el.ref_value * k;
    
    el.unit = simulation_units.current;
    el.si_unit = si_units.current;
    
    circuit.add_element(el);
end

end