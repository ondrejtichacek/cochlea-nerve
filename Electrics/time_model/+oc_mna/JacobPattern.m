function jac_pattern = MNACochleaJacobPattern(A0, channels, Numstacks, ...
    g_BM_displ, I_BM_displ, g_OHC_cilia, I_OHC_cilia, ...
    g_vihc, dg_vihc, I_vihc, g_vohc, dg_vohc, I_vohc, xpos)
% MNACOCHLEAJACOBPATTERN Calculates sparsity pattern of Jacobian
%
% J = df/dy
%
% where C dy/dt = -A(t,y)*y + z = f(x,t)
%
% if A(t,y) = A(t)
% then J = -A(t)
%
% but otherwise, it is a bit more complicated:
% J = -A(t,y) - dA(t,y)/dy * y
%

% Dummy voltage
VoltageIHC = ones(size(xpos));
VoltageOHC = ones(size(xpos));

% Dummy BM displacement
BM_displ = zeros(size(xpos));
OHC_cilia = zeros(size(xpos));

% Modify g
for i = 1:numel(g_BM_displ)
    g_BM_displ{i} = @(y) ones(Numstacks,1);
end

% Modify g
for i = 1:numel(g_OHC_cilia)
    g_OHC_cilia{i} = @(y) ones(Numstacks,1);
end

% Modify g
for i = 1:numel(g_vihc)
    g_vihc{i} = @(V) ones(Numstacks,1);
end

% Modify g
for i = 1:numel(g_vohc)
    g_vohc{i} = @(x,V) ones(Numstacks,1);
end

% Calculate A matrix
A = A_matrix(A0, Numstacks, ...
    g_BM_displ, I_BM_displ, g_OHC_cilia, I_OHC_cilia, ...
    g_vihc, I_vihc, g_vohc, I_vohc, xpos, ...
    BM_displ, OHC_cilia, VoltageIHC, VoltageOHC);


jac_pattern = sparse(A ~= 0);

%% Extension

% this is wrong 

[~, ~, ~, JJ] = ChannelspOpen(struct('IHC', VoltageIHC, 'OHC', VoltageOHC), channels, Numstacks);
[a,b] = size(A);
[aa,bb] = size(JJ);

jac =  [      A, sparse(a,bb);
     sparse(aa,b),          JJ];

 
jac_pattern = sparse(jac ~= 0);

end