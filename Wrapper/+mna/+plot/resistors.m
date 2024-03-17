function [ ] = resistors( xx, VoltIHC_0, VoltOHC_0, mnaopt, plotopt, runopt, unit )
%RESISTORS 

% steady state value of BM and OHC stereocilia displacement
BM_displ = zeros(mnaopt.Numstacks,1);
OHC_cilia = zeros(mnaopt.Numstacks,1);

circuit = mnaopt.circuit;

vi_steady_state = mnaopt.IHC_V_rest(xx);
vo_steady_state = mnaopt.OHC_V_rest(xx);

VoltIHC_0 = VoltIHC_0(:);
VoltOHC_0 = VoltOHC_0(:);

p05 = 0.5*ones(size(vo_steady_state));

deps_1 = { ...
        'xpos', xx, ...
        'BM_displ', BM_displ, ...
        'OHC_cilia', OHC_cilia, ... 
        'vihc', vi_steady_state, ...
        'vohc', vo_steady_state, ...
        'vihc_ss', vi_steady_state, ...
        'vohc_ss', vo_steady_state};

deps_2 = { ...
        'xpos', xx, ...
        'BM_displ', BM_displ, ...
        'OHC_cilia', OHC_cilia, ... 
        'vihc', VoltIHC_0, ...
        'vohc', VoltOHC_0, ...
        'vihc_ss', VoltIHC_0, ...
        'vohc_ss', VoltOHC_0};

channels = mnaopt.IHC_basolateral_channels;

for i = 1:numel(channels)
    p = ['popen_', channels(i).name];

    V = vi_steady_state;
    [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(i));
    deps_1 = [deps_1, {p, Oinf}];

    V = VoltIHC_0;
    [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(i));
    deps_2 = [deps_2, {p, Oinf}];
end

channels = mnaopt.OHC_basolateral_channels;

for i = 1:numel(channels)
    p = ['popen_', channels(i).name];

    V = vo_steady_state;
    [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(i));
    deps_1 = [deps_1, {p, Oinf}];

    V = VoltOHC_0;
    [~, ~, ~, ~, ~, ~, ~, Oinf] = pOpenIHC(V, channels(i));
    deps_2 = [deps_2, {p, Oinf}];
end

selections = { ...
    'type resistor and special radial', ...
    'type resistor and special longitudinal', ...
    'type resistor and special not radial and special not longitudinal', ...
    };

switch unit
    case 'siemens'
        val = @(v) 1e9 * v;
        y_label = 'Conductance (nS)';
    case 'ohm'
        val = @(v) 1 ./ v ./ 1e6;
        y_label = 'Conductance (MOhm)';
    otherwise
        error('not set up')
end
        
for j = 1:numel(selections)
    sel = selections{j};
    ele = circuit.select(sel);

    if isempty(ele)
        continue
    end
    
    assert(all(strcmp(ele(1).unit, {ele.unit})), ...
        'The units of all resistors must be same for plotting');
    
    hfig = figure;
    hold on
    cmp = rwb(numel(ele));
    
    for i = 1:numel(ele)
        ele(i)
        
        assert(strcmp(ele(i).si_unit, 'S'))
        
        % if strcmp(circuit.resistors(1).si_unit, 'S')
        %     ylabel(sprintf('Conductance [n%s]', circuit.resistors(1).unit));
        % elseif strcmp(circuit.resistors(1).si_unit, 'Ohm')
        %     ylabel(sprintf('Resistance [n%s]', circuit.resistors(1).unit));
        % end
        
        v1 = val(ele(i).value(deps_1{:}));
        v2 = val(ele(i).value(deps_2{:}));
        
        plot(xx, v1, 'LineStyle', ':', 'Color', cmp(i,:))
        plot(xx, v2, 'LineStyle', '-', 'Color', cmp(i,:))
        legendtext{2*i-1} = sprintf('Target %s', ele(i).name);
        legendtext{2*i} = sprintf('Value %s', ele(i).name);
    end
    
    ylabel(y_label);
    
    l = legend(legendtext);
    l.ItemHitFcn = @plotOpt.emph_line_width;
    set(gca,'YScale', 'log')

    xlabel('Distance from stapes (norm.)');
   
    if runopt.save_figures
        mySaveFig(hfig, sprintf('resistors_%d_BM_scaling', j), [], runopt, ...
            'relativeDataPath', 'img/MNA/');
    end

end
    
%%

if false
    hfig = figure;
    hold on
    
    RAI = circuit.get_element_by_name('RAI');
    RBI_f = circuit.get_element_by_name('RBI_IHC_fast');
    RBI_s = circuit.get_element_by_name('RBI_IHC_slow');
    
    vA = val(RAI.value(deps_2{:}));
    vB = val(RBI_f.value(deps_2{:})) + val(RBI_s.value(deps_2{:}));
    
    plot(xx, vA ./ vB, ...
        'DisplayName', 'RAI / RBI', ...
        'LineStyle', '-', ...
        'Color', 'b')
    
    RAO = circuit.get_element_by_name('RAO');
    RBO_f = circuit.get_element_by_name('RBO_fast');
    RBO_s = circuit.get_element_by_name('RBO_slow');
    
    vA = val(RAO.value(deps_2{:}));
    vB = val(RBO_f.value(deps_2{:})) + val(RBO_s.value(deps_2{:}));
    
    plot(xx, vA ./ vB, ...
        'DisplayName', 'RAO / RBO', ...
        'LineStyle', '-', ...
        'Color', 'r')
    
    
    ylabel(y_label);
    
    l = legend();
    l.ItemHitFcn = @plotOpt.emph_line_width;
    % set(gca,'YScale', 'log')
    
    xlabel('Distance from stapes (norm.)');
    
    if runopt.save_figures
        mySaveFig(hfig, sprintf('resistors_%d_BM_scaling', j), [], runopt, ...
            'relativeDataPath', 'img/MNA/');
    end
end

end

