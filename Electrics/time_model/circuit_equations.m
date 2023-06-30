function [ C, C0, G, G0, Z, Y0, p, q, vihc, vohc ] = circuit_equations(mnaopt, args)
%CIRCUIT_EQUATIONS calls CIRCUIT_MATRICES
% until convergence is reached between the input and output
% resting potentials, vihc & vohc and Y0.
%
% REMARK 1: the variables p, q are index arrays. They are present so that
%           the optimal order of equations could be returned in them. This
%           feature is currently NOT implemented.
%
% REMARK 2: G is not used by MNAcochlea, and the functions used to create
%           it are never called again. The variable resistors are added to
%           G0 elsewhere.
arguments
    mnaopt (1,1) mnaOpt
    args.initial_guess (1,1) struct = struct();
    args.plotflag (1,1) logical = false
end


% steady state value of BM and OHC stereocilia displacement
BM_displ = zeros(mnaopt.Numstacks,1);
OHC_cilia = zeros(mnaopt.Numstacks,1);

% resting potentials of the inner and outer hair cells before iterative setting
if isempty(fieldnames(args.initial_guess))
    vihc = mnaopt.IHC_V_rest(mnaopt.xgrid(:));
    vohc = mnaopt.OHC_V_rest(mnaopt.xgrid(:));
else
    vihc = args.initial_guess.vihc;
    vohc = args.initial_guess.vohc;
end

% unpack the circuit elements
Resistor = mnaopt.circuit.resistors;
Capacitor = mnaopt.circuit.capacitors;
Vsource = mnaopt.circuit.vsources;
Isource = mnaopt.circuit.isources;

% number of each circuit element
numRes = mnaopt.circuit.num_res;
numCap = mnaopt.circuit.num_cap;
numV = mnaopt.circuit.num_v;
numI = mnaopt.circuit.num_i;

% number of nodes in the circuit
numNode = mnaopt.circuit.num_node;

% some tolerance threshold
volt_EPS = 1e-8 * Unit.conversionConstant('mV', mnaopt.simulation_units.voltage);

store_matrices_as_sparse = true;
if ~store_matrices_as_sparse
    warning('Circuit matrices are not sparse - use this for debug only!')
end

if args.plotflag == true
    figure
    hold on
end

% circuit_matrices() uses vihc and vohc as input to produce Y0, 
% that contains among others the initial nodal voltages.
% The value of vihc and vohc can be recalculated from Y0.
% The loop below performs the modification of vihc and vohc
% until convergence is reached with respect to Y0.

DEBUG = false;
% DEBUG = true;

if DEBUG == true
    figure
    hold on
    plot(vihc*1e3, 'k')
end

Y0_init = 0;

channels = mnaopt.channels;

channels_var_name = {};
for j = 1:numel(channels)
    channels_var_name{j} = ['popen_', channels(j).name];
end

do.mass = false;
do.system = false;
do.rhs = true;
do.jacobian = false;
do = namedargs2cell(do);

for i = 1:100
    
    V = struct('IHC', vihc, 'OHC', vohc, 'IHC_ss', vihc, 'OHC_ss', vohc);
    
    Oinf = {};
    for j = 1:numel(channels)
        [~, ~, ~, ~, ~, ~, ~, Oinf{j}] = pOpenIHC(V.(channels(j).voltage), channels(j), do{:});
        Oinf{j} = num2cell(Oinf{j});
    end
    
    tmp = [channels_var_name; Oinf];
    
    dependencies = struct( ...
        'xpos', num2cell(mnaopt.xgrid(:)), ...
        'BM_displ', num2cell(BM_displ), ...
        'OHC_cilia', num2cell(OHC_cilia), ... 
        'vihc', num2cell(vihc), ...
        'vohc', num2cell(vohc), ...
        'vihc_ss', num2cell(vihc), ...
        'vohc_ss', num2cell(vohc), ...
        tmp{:});
    
    
    [C, C0, G, G0, Z, Y0] = circuit_matrices( ...
        mnaopt.xgrid, ...
        Resistor, Capacitor, Vsource, Isource, ...
        numRes, numCap, numV, numI, numNode, ...
        mnaopt.Numstacks, mnaopt.NumstacksScaling, ...
        dependencies, store_matrices_as_sparse);

    % norm(Y0 - Y0_init);
    Y0_init = Y0;
    
    VI = getVoltage(Y0', 'IHC', mnaopt);
    VO = getVoltage(Y0', 'OHC', mnaopt);
    
    VI = VI(:);
    VO = VO(:);
    
    if args.plotflag == true
        plot(mnaopt.xgrid, VI);
        plot(mnaopt.xgrid, VO);
    end
    
    n1 = norm(VI - vihc) / mnaopt.Numstacks;
    n2 = norm(VO - vohc) / mnaopt.Numstacks;
    
    if n1 < volt_EPS && n2 < volt_EPS
        break
    end

    if DEBUG == true
        plot(VI*1e3)
%         plot(VO*1e3)
    end
    
    vihc = (VI + vihc)/2;% - 0.002;
    vohc = (VO + vohc)/2;
    
%     if DEBUG == true
%         plot(vihc*1e3)
%         plot(vohc*1e3)
%     end
    
end

for j = 1:numel(channels)
    [~, ~, ~, ~, ~, ~, ~, Oinf{j}] = pOpenIHC(V.(channels(j).voltage), channels(j), do{:});
    Oinf{j} = num2cell(Oinf{j});
end

tmp = [channels_var_name; Oinf];

dependencies = struct( ...
    'xpos', num2cell(mnaopt.xgrid(:)), ...
    'BM_displ', num2cell(BM_displ), ...
    'OHC_cilia', num2cell(OHC_cilia), ... 
    'vihc', num2cell(vihc), ...
    'vohc', num2cell(vohc), ...
    'vihc_ss', num2cell(vihc), ...
    'vohc_ss', num2cell(vohc), ...
    tmp{:});


[C, C0, G, G0, Z, Y0] = circuit_matrices( ...
    mnaopt.xgrid, ...
    Resistor, Capacitor, Vsource, Isource, ...
    numRes, numCap, numV, numI, numNode, ...
    mnaopt.Numstacks, mnaopt.NumstacksScaling, ...
    dependencies, store_matrices_as_sparse);

% VI = getVoltage(Y0', 'IHC', mnaopt);

% permutation is not yet fully implemented
% remaining:
%    - y tolerance for ode15s conditional of y index
%    - Jacobian - depends on y
%    - ?

% [p,q,r,s,cc,rr] = dmperm(G);

% [p,q,r] = btf(G);


[p, q] = deal(1:numel(Z));

% p = randperm(numel(Z));
% q = randperm(numel(Z));

% hack sparse allocation of C0 and A0
[a,b] = size(C0);

tmp = spalloc(a,b,nnz(C));
tmp(1:a,1:b) = C0;
C0 = tmp;

assert(nzmax(C0) == nnz(C));

[a,b] = size(G0);

tmp = spalloc(a,b,nnz(G));
tmp(1:a,1:b) = G0(1:a,1:b);
G0 = tmp;

assert(nzmax(G0) == nnz(G));

% C0 = C0(p,q); % do not do
C = C(p,q);
G = G(p,q);
% G0 = G0(p,q); % do not do
Y0 = Y0(q);
Z = Z(p);

% i = [9,2,3,4,5,6,7,8,1];
% 
% C_ = C(i,:);
% G_ = G(i,:);
% G0_ = G0(i,:);
% Z_ = Z(i);
% Y0_ = Y0(i);
% 
% n = numel(Z);
% 
% N = 100000;
% r = zeros(1,N);
% rbest = 0;
% 
% for i = 1:N
%     I = randperm(n);
%     r(i) = 1/condest(G(I,:),10);
%     if r(i) > rbest
%         rbest = r(i);
%         Iberst = I;
%     end
% end

end