function [ y_out, t_out, V_out, y_dyn_out, n_out, tropt_out, f_ChannelsOpenSS_out, f_CalciumCurrent_out, f_TransmitterRelease_out ] = Transduction_v4_multi( ...
        t, IHCVoltage, SR, replica, topt, antopt, runopt, memopt )
%TRANSDUCTION

x_pos = antopt.positionsToExcite;

[y_out, t_out, V_out, y_dyn_out, n_out, tropt_out, f_ChannelsOpenSS_out, f_CalciumCurrent_out, f_TransmitterRelease_out] = ...
    deal(cell(numel(x_pos), 1));

for i = 1:numel(x_pos)
    antopt_cpy = copy(antopt);
    antopt_cpy.slicesToExcite = i;
    antopt_cpy.positionsToExcite = x_pos(i);
    V = IHCVoltage(:,i);
    [y_file, t_file, V_file, y_dyn_file, n_out{i}, tropt_out{i}, f_ChannelsOpenSS_out{i}, f_CalciumCurrent_out{i}, f_TransmitterRelease_out{i} ] = ...
        Transduction_v4( t, V, SR, replica, topt, antopt_cpy, runopt, memopt );

    yy = load(y_file.Properties.Source);
    tt = load(t_file.Properties.Source);
    VV = load(V_file.Properties.Source);
    yydd = load(y_dyn_file.Properties.Source);

    y_out{i} = yy;
    t_out{i} = tt;
    V_out{i} = VV;
    y_dyn_out{i} = yydd;
end

tmp = [y_out{:}];
y_file.y = cat(3, tmp.y);
y_file.vesicle_release_events = {tmp.vesicle_release_events};


tmp = [t_out{:}];
t_file.t = cat(2, tmp.t);
tmp = [V_out{:}];
V_file.V = cat(2, tmp.V);

tmp = [y_dyn_out{:}];
y_dyn_file.y = cat(2, tmp.y);

y_out = y_file;
t_out = t_file;
V_out = V_file;
y_dyn_out = y_dyn_file;


end