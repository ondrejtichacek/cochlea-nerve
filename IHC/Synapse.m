function [ synapseResult ] = Synapse( V, T, antopt, sr, replica, topt, runopt, memopt )
%SYNAPSE Summary of this function goes here
%
%   Voltage = IHC voltage
%   Time    = time points
%   hhopt   = options for Hodgkin-Huxley
%   sr      = spontaneous rate (name)

mkdir(replica.dir);

if isempty(replica.rng)
    rng(replica.seed, 'twister');
else
    rng(replica.rng)
end

[Time, ind] = replica.tspan_complement.intersect_time(T.Value);

Voltage = V.Value(ind,:);

switch antopt.ant

    case {'v4'}
        Voltage = Voltage *1e-3; % [mV] to [V]

        [y, t, V, y_dyn, n, tropt, f_ChannelsOpenSS, f_CalciumCurrent, f_TransmitterRelease ] = ...
            Transduction_v4_multi( Time, Voltage, sr, replica, topt, antopt, runopt, memopt );

    otherwise
        error('ANT:InvalidArguments', ['Unsupported option for parametr hhopt.ant: ', antopt.ant]);
        
end

%% store results

rng_state = rng;

synapseResult = SynapseResult();

synapseResult.protect_PATH_error = replica.protect_PATH_error;

synapseResult.PATH = replica.dir;
synapseResult.orig_PATH = replica.dir;

synapseResult.tropt = tropt;

synapseResult.fy = y;                 % ... MatFile
synapseResult.ft = t;                 % ... MatFile
synapseResult.fV = V;                 % ... MatFile
synapseResult.fy_dyn = y_dyn;          % ... MatFile
synapseResult.ChannelsOpenSS     = f_ChannelsOpenSS;     % ... function handle
synapseResult.CalciumCurrent     = f_CalciumCurrent;     % ... function handle
synapseResult.TransmitterRelease = f_TransmitterRelease; % ... function handle


synapseResult.size_info_dyn = y_dyn.size_info;
synapseResult.size_info = y.size_info;

synapseResult.name = antopt.ant;    % ... string
synapseResult.n = n;                  % ... double
synapseResult.n_dyn_rep = numel(y_dyn.y);

synapseResult.xgrid_ind = antopt.slicesToExcite;
synapseResult.xgrid = antopt.positionsToExcite;

synapseResult.SR = sr;  % ... string

synapseResult.rng_state = rng_state;
synapseResult.seed = rng_state.Seed;

save(fullfile(replica.dir, 'synapseResult.mat'), 'synapseResult', '-v7.3')


end

