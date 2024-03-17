function [mnaopt] = ocelectric_templates(mnaopt, args)
arguments
    mnaopt (1,1) mnaOpt
    args.mna_ver (1,:) char = 'dev2'
end

fprintf("*** Using %s version of MNA ***\n", args.mna_ver);

kwargs = {};

switch args.mna_ver    
    case 'old'
        mnaopt.capacitance_state_dependence = 'none';
        mnaopt.IHC_basolateral_conductance_dependance = 'vihc_ss';
        mnaopt.OHC_basolateral_conductance_dependance = 'vohc_ss';
        OHC_NLC_description = 'linear';
        circuit_connection = 'v2';
    case 'new'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        mnaopt.OHC_basolateral_conductance_dependance = 'vohc';
        OHC_NLC_description = 'nonlinear_freq';
    case 'dev'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3';
    case 'dev2_ss_nlc'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq_ss';
        circuit_connection = 'v3';
    case 'dev2_ss_nlc_v2'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq_ss_v2';
        circuit_connection = 'v3';
    case 'dev2_no_nlc'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'linear';
        circuit_connection = 'v3';
    case 'dev2'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3';
    case 'dev2_ss'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        mnaopt.IHC_basolateral_popen_dependance = 'vihc_ss';
        mnaopt.OHC_basolateral_popen_dependance = 'vohc_ss';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3';
    case 'dev2_ss_ihc'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        mnaopt.IHC_basolateral_popen_dependance = 'vihc_ss';
        % mnaopt.OHC_basolateral_popen_dependance = 'vohc_ss';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3';
    case 'dev2_ss_ohc'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        % mnaopt.IHC_basolateral_popen_dependance = 'vihc_ss';
        mnaopt.OHC_basolateral_popen_dependance = 'vohc_ss';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3';
    case 'dev2_IHC_stereocilia_noise_damage'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3';
        kwargs = {
            'IHC_stereocilia_noise_damage', true, ...
            'IHC_stereocilia_noise_damage_f0', Frequency(2400, 'Hz'), ...
            'IHC_stereocilia_noise_damage_degree', 0.01, ...
            };
    case 'dev2_OHC_stereocilia_noise_damage'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3';
        kwargs = {
            'OHC_stereocilia_noise_damage', true, ...
            'OHC_stereocilia_noise_damage_f0', Frequency(2400, 'Hz'), ...
            'OHC_stereocilia_noise_damage_degree', 0.01, ...
            };
    case 'dev3'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v3r';
    case 'dev_v0e'
        mnaopt.capacitance_state_dependence = 'strong';
        mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_tonotopic';
        % mnaopt.IHC_basolateral_conductance_dependance = 'vihc';
        % mnaopt.IHC_basolateral_conductance_dependance = 'channel_popen_Dierich_2020';
        mnaopt.OHC_basolateral_conductance_dependance = 'channel_popen';
        OHC_NLC_description = 'nonlinear_freq';
        circuit_connection = 'v0e';
end

mnaopt.circuit = cochlear_circuit(mnaopt, false, ...
    ...'OHC_NLC_description', 'linear', ...
    ...'OHC_NLC_description', 'nonlinear', ...
    ...'OHC_NLC_description', 'nonlinear_freq', ...
    'OHC_NLC_description', OHC_NLC_description, ...
    ...'connection', 'mistrik');
    'connection', circuit_connection, ...
    kwargs{:});

end