function [ g, dg, Ig, c, Ic, dep_fun, dep_sel, dep_sign] ...
        = JacobPrecalc( Circuit, NumstacksScaling, x, mnaopt )
arguments
    Circuit
    NumstacksScaling
    x
    mnaopt (1,1) mnaOpt
end
% JACOBPRECALC Pre-calculates (circuit related) constants used in the MNACochleaJacob function

channels = mnaopt.channels;

Resistor = Circuit.resistors;
Capacitor = Circuit.capacitors;

[nn, mm] = deal(Circuit.num_v + Circuit.num_node);

numNode = Circuit.num_node;

[g, dg, Ig, c, Ic] = deal(struct());

deps = {'BM_displ', 'OHC_cilia', 'vihc', 'vohc', 'vihc_ss', 'vohc_ss'};
if ~isempty(channels)
    deps = [deps, channels.var_name];
end

for dep = string(deps)
    [g.(dep), dg.(dep), Ig.(dep)] = conductance_function( ...
            Resistor, dep, numNode, nn, mm, x, NumstacksScaling);
    
    [c.(dep), ~, Ic.(dep)] = conductance_function( ...
            Capacitor, dep, numNode, nn, mm, x, NumstacksScaling);
end

%% IHC & OHC voltage

[dep_fun.vihc, dep_sel.vihc, dep_sign.vihc] = get_voltage_from_selection('IHC', mnaopt);
[dep_fun.vohc, dep_sel.vohc, dep_sign.vohc] = get_voltage_from_selection('OHC', mnaopt);

for k = 1:numel(channels)
    channel = channels(k);
    [~, ~, dep_fun.(channel.var_name)] = getChannel([], channel, mnaopt);
end

%%

    function [sel_function, sel, sign_part] = get_voltage_from_selection(name, mnaopt)
        [~, sel, ~, operand, ~] = getVoltage([], name, mnaopt);

%         assert(operand == '-');
        assert(numel(sel) == 2);
        
        switch operand
            case '-'
                sel_function = @(y) y(sel{1}) - y(sel{2});
                sign_part = [+1, -1];
            case '+'
                sel_function = @(y) y(sel{1}) + y(sel{2});
                sign_part = [+1, +1];
        end
    end


    function [g, dg, I] = conductance_function(Element, dependence, numNode, nn, mm, x, NumstacksScaling)
        % extract BM-dependent resistors
        E = Element(Element.depends_on(dependence));

        % prealoc
        [I,g,dg] = deal(cell(1,numel(E)));        

        % compute
        for i = 1:numel(E)
            
            E(i).debug.value_assigned = dependence;
            
            if isempty(E(i).pre)
                r = E(i).value_fun('xpos', x);
                [dr.fun, dr.args] = E(i).der_fun({dependence}, 'xpos', x);
            else
                r = E(i).value_fun_pre;
                if isfield(E(i).pre, 'dep_fun_der')
                    [dr.fun, dr.args] = E(i).der_fun_pre;
                else
                    [dr.fun, dr.args] = E(i).der_fun({dependence}, 'xpos', x);
                end
            end
            
            switch Element(i).type
                case 'resistor'
                    switch Element(i).si_unit
                        case 'S'
                            g{i} = r;
                            dg{i} = dr;
                        case 'Ohm'
                            g{i} = @(y) 1./ r(y);
                            dg{i} = @(y) 1./ dr(y);
                        otherwise
                            error('Unknown unit %s.', Element(i).si_unit);
                    end
                case 'capacitor'
                    g{i} = r;
                    dg{i} = dr;
                otherwise
                    error('Unsupported element type %s.', Element(i).type);                    
            end
            
            % G(t) PATTERN
            % I containts a pattern where time-dependent values of circuit are stored
            % in the matrix. It contains a sign, but not value. Value is stored in the
            % vector g. The third dimmension of I corresponds to the number of time
            % dependent elements.

            assert(any(strcmp(Element(i).type, {'resistor','capacitor'})))
            I{i} = sparse(NumstacksScaling * Gt_pattern(E(i), numNode, nn, mm));
            
            v = zeros(1,size(I{i},1));
            for j = 1:size(I{i},1)
                v(j) = any(I{i}(j,:));
            end   
            v = find(v);
            
            dg{i}.nodes_affected = v;
            
        end
    end
end