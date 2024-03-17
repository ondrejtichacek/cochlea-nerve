function [ hfig ] = voltage_mesh_one( Xgrid, Tgrid, Volt, t0, tf, plotopt, V0)
arguments
    Xgrid (:,:) double
    Tgrid (:,:) double
    Volt (:,:) double
    t0 (1,1) Time
    tf (1,1) Time
    plotopt
    V0 = []
end
%VOLTAGE_MESH_ONE

hfig = plotopt.figure;

if isempty(V0)
    V0 = Volt(1,:);
end

C = Volt - V0;

M = max(abs(C(:)));

surf(Xgrid, Tgrid, Volt, C, ...
    plotopt.surf_options{:});

xlabel('BM position');
ylabel('Time (ms)');
zlabel('Voltage (mV)');

if ~isempty(t0) && ~isempty(tf)
    ylim([t0.ms, tf.ms]);
end

% M = max(abs(Volt(:) - V0));
% caxis([V0-M, V0+M])

caxis([-M, M])

colormap(plotopt.centered_colormap())
view(plotopt.view);

set(gca, 'color', 0.95*ones(1,3))

clear Volt

end

