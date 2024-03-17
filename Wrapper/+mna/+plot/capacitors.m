function [ ] = capacitors( xx, VoltIHC_0, VoltOHC_0, mnaopt, plotopt, runopt )
%CAPACITORS 

% steady state value of BM and OHC stereocilia displacement
BM_displ = zeros(mnaopt.Numstacks,1);
OHC_cilia = zeros(mnaopt.Numstacks,1);

circuit = mnaopt.circuit;

vi_steady_state = mnaopt.IHC_V_rest(xx);
vo_steady_state = mnaopt.OHC_V_rest(xx);

VoltIHC_0 = VoltIHC_0(:);
VoltOHC_0 = VoltOHC_0(:);

deps = { ...
        'xpos', xx, ...
        'BM_displ', BM_displ, ...
        'OHC_cilia', OHC_cilia, ... 
        'vihc', VoltIHC_0, ...
        'vohc', VoltOHC_0, ...
        'vihc_ss', VoltIHC_0, ...
        'vohc_ss', VoltOHC_0};


hfig = figure;
    

hold on
for i = 1:numel(circuit.capacitors)
    plot(xx,circuit.capacitors(i).value(deps{:}))
    legendtext{i} = circuit.capacitors(i).name;
end
l = legend(legendtext);
l.ItemHitFcn = @plotOpt.emph_line_width;
set(gca,'YScale', 'log')

xlabel('Distance from stapes (norm.)');
ylabel('Capacitance [F]');

if runopt.save_figures
    mySaveFig(hfig, 'capacitors_BM_scaling', [], runopt, ...
        'width', '0.45\textwidth', ...
        'height', '0.5\textwidth', ...
        'relativeDataPath', 'img/MNA/');
end

end

