function [ f, B, z ] = System( t, y, ...
        z, B0, channels, Numstacks, size_info, mechfs, ...
        G, dG, I, xpos, dep_fun, dep_sel, dep_sign, ...
        dep_names, dep_names_alt, ...
        p, q, needs_ordering, ...
        mechopt, t_mech, loaded_BM_displ, loaded_OHC_cilia)
arguments
    t (1,1) double
    y (:,1) double
    z (:,1) double
    B0 (:,:) double
    channels (1,:) struct
    Numstacks (1,1) double
    size_info (1,1) struct
    mechfs (1,1) double
    G (1,1) struct
    dG (1,1) struct
    I (1,1) struct
    xpos (:,1) double
    dep_fun (1,1) struct
    dep_sel (1,1) struct
    dep_sign (1,1) struct
    dep_names (1,:) cell
    dep_names_alt (1,:) cell
    p (1,:) double
    q (1,:) double
    needs_ordering (1,1) logical
    mechopt (1,1) mechOpt
    t_mech double = []
    loaded_BM_displ double = []
    loaded_OHC_cilia double = []
end
%DERIVATION: Calculates RHS of equation
%
%  C dy/dt = -A(t,y)*y + z(t,y) = f(x,t)
%          =: B(t,y)*y + z(t,y)

if needs_ordering
    % inverse permutations
    P(p) = 1:numel(p);
    Q(q) = 1:numel(q);

    yQ = y(Q);
else
    yQ = y;
end

% Interpolate BM displacement
switch mechopt.integration
    case 'standalone'
        [BM_displ, OHC_cilia] = BMdispl(t, t_mech, mechfs, loaded_BM_displ, loaded_OHC_cilia);

    case 'electric'
        [BM_displ, OHC_cilia] = mech_interpolate(t, y, xpos, size_info, mechopt);
    otherwise
        error('Unsupported value of mechopt.integration = %s', mechopt.integration)
end

vars = struct( ...
    'BM_displ', BM_displ(:), ...
    'OHC_cilia', OHC_cilia(:) ...
    );

% Extract IHC & OHC receptor potentials from solution
fn = fieldnames(dep_fun);
for i = 1:numel(fn)
    dep = fn{i};
    vars.(dep) = dep_fun.(dep)(yQ);
end

% Link to steady-state solution
vars.vihc_ss = mechopt.vihc_ss;
vars.vohc_ss = mechopt.vohc_ss;

% Calculate A matrix from the constant A0
B = oc_mna.SystemMatrix(B0, Numstacks, G, I, xpos, vars, size_info);

%%
if ~isempty(channels)

    V = struct('IHC', vars.vihc, 'OHC', vars.vohc, 'IHC_ss', vars.vihc_ss, 'OHC_ss', vars.vohc_ss);

    [~, AA, zz] = ChannelspOpen(V, channels, Numstacks, ...
        'mass', false, 'system', true, 'rhs', true, 'jacobian', false);

    % [a,b] = size(A);
    % [aa,bb] = size(AA);
    % A =  [         A, sparse(a,bb);
    %      sparse(aa,b),         AA];

    m = size_info.channels.start;
    n = size_info.channels.end;

    B(m:n,m:n) = -AA;

    z = [z; zz];
end

%% Extend system by F = FOHC ~ high pass filtered VOHC

if true%false
    z = [z; zeros(mechopt.Numstacks, 1)];
end

%% Mechanical PDE
% M * dy/dt = J * y + g(t) + h(y)
if strcmp(mechopt.integration, 'electric')

    n = numel(z) + 1;

    % A_mech is already contained in A0 since it's constant
    % A_mech = - mechopt.jacobian;
    y_mech = y(n:end);

    g = mechopt.OdeFun_t;

    switch mechopt.amplifier
        case 'none'

            z_mech = g(t);

        case 'mechanic'

            h = mechopt.OdeFun_y;
            z_mech = g(t) + h(y_mech);

        case 'electric'

            if isnan(mechopt.vohc_0)
                eohc = zeros(size(vars.vohc));
            else
                eohc = vars.vohc - mechopt.vohc_0;
            end

            h = mechopt.OdeFun_v;

            z_mech = g(t) + h(eohc);
    end

    % z_mech = - z_mech;

    % [a,b] = size(A);
    % [aa,bb] = size(A_mech);
    % 
    % A =  [         A, sparse(a,bb);
    %      sparse(aa,b),    A_mech];

    z = [z; z_mech];
end

%% Extend system by regularization term

if true%false
    % reg = 0;
    reg = V.OHC(1) > -20e-3 || V.OHC(1) < -120e-3;
    z = [z; reg];
end

%%

% figure(111)
% hold on
% plot(y(size_info.channel_OHC_fast.start:2:size_info.channel_OHC_fast.end))
% plot(y(size_info.channel_OHC_slow.start:2:size_info.channel_OHC_slow.end))


% ordering
if needs_ordering
    z = z(q);
    B = B(p,q);
end

% assemble f(x,t)
[n, m] = size(y);

if m > 1 % prealocation
    f = zeros(n,m);
end

for i = 1:m
    f(:,i) = B * y(:,i) + z;
end

end