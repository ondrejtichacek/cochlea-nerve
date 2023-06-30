function [ AA, CC, ZZ, Y_desc ] = RLC( ...
        Resistor, Inductor, Capacitor, Vsource, Isource, numNode)

% [ C 0 0 ]      [ v   ]   [  G     A_L  A_V ] [v  ]   [ -A_I * I ]
% [ 0 L 0 ] d/dt [ i_L ] + [ -A_L'  0    0   ] [i_L] = [ 0        ]
% [ 0 0 0 ]      [ i_V ]   [  A_V'  0    0   ] [i_V]   [ E        ]

%  CC       d/dt x       + AA x                      = ZZ 

    
%%

numR = numel(Resistor);
numL = numel(Inductor);
numC = numel(Capacitor);
numV = numel(Vsource);
numI = numel(Isource);


%% Preallocate

G = zeros(numNode, numNode);    % G is the resistor matrix
F = zeros(numNode, numNode);    % F is the capacitor matrix

BV = zeros(numNode, numV);
CV = zeros(numV, numNode);
DVV = zeros(numV, numV);

BL = zeros(numNode, numL);
CL = zeros(numL, numNode);
DLL = zeros(numL, numL);

DLV = zeros(numL, numV);
DVL = zeros(numV, numL);

I = zeros(numNode, 1);       % current sources
E = zeros(numV, 1);          % voltage sources
H = zeros(numL, 1);

KV = zeros(numV, numNode);
JV = zeros(numNode, numV);

KL = zeros(numL, numNode);
JL = zeros(numNode, numL);

L = zeros(numL, numL);         % L is the inductor matrix

NVL = zeros(numV, numL);
NLV = zeros(numL, numV);
NVV = zeros(numV, numV);


%% Fill the G matrix with resistros

for i = 1:numR
    n1 = Resistor(i).node1;
    n2 = Resistor(i).node2;
    
    g = Resistor(i).value();
    
    % If neither side of the element is connected to ground
    % then subtract it from appropriate location in matrix.
    if (n1 ~= 0) && (n2 ~= 0)
        G(n1,n2) = G(n1,n2) - g;
        G(n2,n1) = G(n2,n1) - g;
    end
    
    % If node 1 is connected to ground, add element to diagonal of matrix.
    if n1 ~= 0
        G(n1,n1) = G(n1,n1) + g;
    end
    % Ditto for node 2.
    if n2 ~= 0
        G(n2,n2) = G(n2,n2) + g;
    end
    
end
%% Fill the F matrix with capacitors

 for i = 1:numC
    n1 = Capacitor(i).node1;
    n2 = Capacitor(i).node2;
    
    f = Capacitor(i).value();
        
    % If neither side of the element is connected to ground
    % then subtract it from appropriate location in matrix.
    if (n1 ~= 0) && (n2 ~= 0)
        F(n1,n2) = F(n1,n2) - f;
        F(n2,n1) = F(n2,n1) - f;
    end
    
    %If node 1 is connected to ground, add element to diagonal of matrix.
    if n1 ~= 0
        F(n1,n1) = F(n1,n1) + f;
    end
    %Ditto for node 2.
    if n2 ~= 0
        F(n2,n2) = F(n2,n2) + f;
    end
    
 end

 %% Fill the L matrix with inductors

 for i = 1:numL
    L(i,i) = Inductor(i).value();
 end

 %% Fill the I matrix
for j = 1:numNode
    for i = 1:numI
        if Isource(i).node1 == j
            I(j) = I(j) - Isource(i).value();
        elseif Isource(i).node2 == j
            I(j) = I(j) + Isource(i).value();
        end
    end
end
    
%% Fill the B matrix
for i = 1:numV
    for j = 1:numNode
        if Vsource(i).node1 == j
            BV(j,i) = 1;
        elseif Vsource(i).node2 == j
            BV(j,i) = -1;
        end
    end
end

for i = 1:numL
    for j = 1:numNode
        if Inductor(i).node1 == j
            BL(j,i) = 1;
        elseif Inductor(i).node2 == j
            BL(j,i) = -1;
        end
    end
end

%% Fill the C matrix
for i = 1:numV
    for j = 1:numNode
        if Vsource(i).node1 == j
            CV(i,j) = 1;
        elseif Vsource(i).node2 == j
            CV(i,j) = -1;
        end
    end
end

for i = 1:numL
    for j = 1:numNode
        if Inductor(i).node1 == j
            CL(i,j) = -1;
        elseif Inductor(i).node2 == j
            CL(i,j) = 1;
        end
    end
end

%% Fill the D matrix
%The D matrix is non-zero only for CCVS and VCVS (not included
%in this simple implementation of SPICE)


%% Fill the E matrix
for i=1:numV
    E(i) = Vsource(i).value();
end
    
%% Form the AA, CC, and ZZ matrices

AA = [G   BL   BV;
      CL  DLL  DLV;
      CV  DVL  DVV];
  
CC = [F   JL   JV;
      KL  L    NLV;
      KV  NVL  NVV];

ZZ = [I;
      H;
      E];

I_desc = {};
for i = 1:size(I,1)
    I_desc{i} = sprintf('V_N%d',i);
end
H_desc = {};
for i = 1:size(H,1)
    H_desc{i} = sprintf('I_L%d',i);
end
E_desc = {};
for i = 1:size(E,1)
    E_desc{i} = sprintf('I_V%d',i);
end
Y_desc = [I_desc, H_desc, E_desc];

end
