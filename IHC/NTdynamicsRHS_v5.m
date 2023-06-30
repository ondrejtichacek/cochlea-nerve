function [ dz ] = NTdynamicsRHS_v5( t, z, ...
    size_info, channels, vesicles, ps, y, l, x, r, G_Ca, C_Ca_background, ...
    channels_open_ss_parameters_normal, channels_open_ss_parameters_burst, ...
    transmitter_release_parameters, ...
    V_steady_state, ...
    dt, C_vesicles)
arguments
    t
    z
    size_info
    channels
    vesicles
    ps
    y
    l
    x
    r
    G_Ca
    C_Ca_background
    channels_open_ss_parameters_normal
    channels_open_ss_parameters_burst
    transmitter_release_parameters
    V_steady_state
    dt
    C_vesicles
end
%TRANSDUCTIONRHS

C_vesicles = C_vesicles(:);

q = decompose_z(z, 'NT_free', size_info);
c = decompose_z(z, 'NT_cleft', size_info);
w = decompose_z(z, 'NT_reprocessing', size_info);
c_proton = decompose_z(z, 'proton_cleft', size_info);

% iteration index
it = ceil(t/dt);

C_vesicles = C_vesicles + C_Ca_background;

%% Transmitter Release and Recycling

[dq, dc, dw, dc_proton, vesicles] = NTdynamicsRHS_v5_core( t, ...
    q, c, w, c_proton, ...
    vesicles, y, l, x, r, ...
    transmitter_release_parameters, ...
    dt, C_vesicles);


%% Build dz

dz = [dq; dc; dw; dc_proton];

if any(isnan(dz))
    error('NaN encountered')
end

end

function v = decompose_z(z, variable, size_info)

si = size_info.(variable);
v = z(si.start : si.end);

end
