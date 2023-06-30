function [ dq, dc, dw, dc_proton, vesicles] = NTdynamicsRHS_v5_core( t, ...
    q, c, w, c_proton, ...
    vesicles, y, l, x, r, ...
    transmitter_release_parameters, ...
    dt, C_vesicles)
arguments
    t
    q
    c
    w
    c_proton
    vesicles
    y
    l
    x
    r
    transmitter_release_parameters
    dt
    C_vesicles
end
%TRANSDUCTIONRHS

C_vesicles = C_vesicles(:);

M = vesicles.num;

%% Transmitter Release and Recycling

kmax = transmitter_release_parameters{1};
% n_Hill = transmitter_release_parameters{2};
% KA_Hill = transmitter_release_parameters{3};

assert(kmax*dt <= 1)

% k1 ... single spot release rate
% k1 = kmax .* Hill_Langmuir_A(C_vesicles, n_Hill, KA_Hill);
k1 = RibbonSynapse_v4.TransmitterRelease(C_vesicles, transmitter_release_parameters{:});

q_old = q;

% released
Nqk = NTTransport(q, k1*dt, M);
q = q - Nqk;

% reprocessed
Nwx = NTReprocessing(~q, floor(w)*x*dt, M);
q = q + Nwx;

% manufactured
NMqy = NTManufacture(~q, y*dt, M);
q = q + NMqy;

% dq = (Nwx + NMqy - Nqk)/dt;
dq = (q - q_old)/dt;
dc = sum(Nqk)/dt - l.*c - r.*c;
dw = r.*c - sum(Nwx)/dt;

dc_proton = sum(Nqk)/dt - 0.9*l.*c_proton;

for i = 1:sum(Nqk)
    vesicles.logReleaseEvent(t);
end

%% Stochastic NT Transport

function N = NTReprocessing(q, rho, n)
    n_empty = sum(q);
    if n_empty > 0 && rho > 0
        N = q & (rand(n,1) < (rho / n_empty));
    else
        N = zeros(size(q));
    end
end
function N = NTManufacture(q, rho, n)
    n_empty = sum(q);
    if n_empty > 0
        N = q & (rand(n,1) < (rho / n_empty));
    else
        N = zeros(size(q));
    end
end
function N = NTTransport(q, rho, n)
    n_occupied = sum(q);
    if n_occupied > 0
        N = q & (rand(n,1) < rho);
    else
        N = zeros(size(q));
    end
end


end
