function [ nerveResult ] = cochleaRunBruce( topt, hhopt, scaleFactor, slices, replica, runopt, memopt )
% COCHLEARUNBRUCE Runs stochastic_HH_fancyinput for each slice of the cochela
%

%% init

synapseResult = replica.src_data;

mkdir(replica.dir)

if isempty(replica.rng)
    rng(replica.seed, 'twister');
else
    rng(replica.rng)
end


%% preprocessing

[T, ind] = replica.tspan_complement.intersect_time(synapseResult.t);

% we assume that if c = 1, then all receptors are bound and the current is
% maximal =>
% I = min(1,c) * I_max;

%% saturation

% Y = ((1./(1 + c1*exp(-y./y1) + c2*exp(-y./y2))) - b);

c1 = 0.7293;
c2 = 1.4974;
% b = 0.30991;

b = 1/(1 + c1 + c2); % so that Y(0) = 0

% considering y1 < y2 and y1 > 0, y2 > 0:
% Y(-infinity) = -b
% Y(infinity) = 1-b

N = 3; % saturates at this point

q = N / (1-b);

y1 = q*0.1139;
y2 = q*0.3736;

Ynl = @(y) q * ((1./(1 + c1*exp(-y./y1) + c2*exp(-y./y2))) - b);

%%

c = synapseResult.c(ind, ':');
c = Ynl(c);

c_arr = synapseResult.c_dyn(ind, ':');
N_arr = numel(c_arr);

c = zeros(size(c,1), size(c,2), N_arr);
for j = 1:N_arr
    c(:,:,j) = Ynl(c_arr{j});    
end

% CurrIHCSynapse = min(1, c) * scaleFactor;
CurrIHCSynapse = c * scaleFactor;

% scaleFactor is in units [uA/c] where [c] is tje unit of normalized 
% NT concentration in the cleft

%%
dt = mean(diff(T));

assert(hhopt.samplingFrequency >= Frequency(100, 'kHz'))


if abs(hhopt.samplingFrequency - 1/dt) > Frequency(1, 'Hz')

    % assert(numel(slices) == 1) % otherwise not implememnted yet

    F = griddedInterpolant(T.us, CurrIHCSynapse);

    t0 = T(1);
    tf = T(end);

    signalLength = tf - t0; % [s]
    numberOfPoints = 1 + ceil(hhopt.samplingFrequency * signalLength);
    T = linspace( t0, tf, numberOfPoints ); % time points
    T = T';

    
    CurrIHCSynapse = F(T.us);
end

if numel(slices) > 1
    error("fix this for narr > 1")
    assert(size(CurrIHCSynapse, 2) == 1)
    sz = size(CurrIHCSynapse);
    sz(2) = [];
    CurrIHCSynapse = reshape(CurrIHCSynapse, sz);
end

%% computation

% this is a quick workaround
tmp_y = cell(N_arr, 1);
tmp_i = cell(N_arr, 1);

for j = 1:N_arr

    if j < N_arr
        hhopt.save_method = 'struct';
    else
        hhopt.save_method = 'matlab_matfile';
    end

    [y, t, i, n, Na_max, calcN_Na, mh, mhinit] = stochHH( ...
        hhopt.method, T, CurrIHCSynapse(:,slices,j), ...
        replica, topt, hhopt, runopt, memopt);
    
    if any(isnan(y.y))
        error('NaN found in result');
    end

    tmp_y{j} = y.y;
    tmp_i{j} = i.Istim;

end

y.y = tmp_y;
i.Istim = tmp_i;

%% store results

rng_state = rng;

nerveResult = NerveResult();

nerveResult.protect_PATH_error = replica.protect_PATH_error;

nerveResult.PATH = replica.dir;
nerveResult.orig_PATH = replica.dir;

nerveResult.n_dyn_rep = N_arr;

nerveResult.fy = y;                 % ... MatFile
nerveResult.ft = t;                 % ... MatFile
nerveResult.fi = i;                 % ... MatFile
nerveResult.calcN_Na = calcN_Na;    % ... function handle
nerveResult.mh = mh;                % ... function handle
nerveResult.mhinit = mhinit;        % ... function handle

nerveResult.name = hhopt.method;    % ... string
nerveResult.n = n;                  % ... double
nerveResult.Na_max = Na_max;        % ... double

nerveResult.SR = synapseResult.SR;  % ... string

nerveResult.xgrid_ind = synapseResult.xgrid_ind;
nerveResult.xgrid = synapseResult.xgrid;

nerveResult.rng_state = rng_state;
nerveResult.seed = rng_state.Seed;

save(fullfile(replica.dir, 'nerveResult.mat'), 'nerveResult', '-v7.3')

end
