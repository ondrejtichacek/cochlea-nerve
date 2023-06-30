function [ Acell, Gcapcell, Zcell ] = RCconductances( ...
        Resistor, Capacitor, Vsource, Isource, ...
        numRes, numCap, numV, numI, numNode, ...
        omit_state_dep_elements, NumstacksScaling, ...
        dependencies)
%RCCONDUCTANCES Construction of conductive matrices
%
%   Acell:      a complete conductance matrix from resistors
%   Gcapcell:   a complete conductance matrix from capacitors
%   Zcell:      a complete matrix of independent sources


%% Preallocate all of the cell arrays
G = zeros(numNode,numNode);    % G is the resistor matrix      (A)
F = zeros(numNode,numNode);    % F is the capacitor matrix     (K)
I = zeros(numNode,1);          % I current sources             (Z) -- stays zero

if numV ~= 0
    B = zeros(numNode,numV);    % conection of voltage sources  (A)
    C = zeros(numV,numNode);    % conection of voltage sources  (A)
    D = zeros(numV,numV);       % conection of voltage sources  (A) -- stays zero
    
    E = zeros(numV,1);          % independent voltage sources   (z)
    K = zeros(numV,numNode);    % conection of voltage sources  (K) -- stays zero
    J = zeros(numNode,numV);    % conection of voltage sources  (K) -- stays zero
    N = zeros(numV,numV);       % conection of voltage sources  (K) -- stazs zero
end

%% Fill the G matrix with transducer resistros from netlist
if omit_state_dep_elements == true
    % omit RAO and RAI elements in building of the G matrix
    G = G0_matrix(Resistor, numRes, numNode, NumstacksScaling, dependencies);
else
    % include all elements as normal
    G = G_matrix(Resistor, numRes, numNode, NumstacksScaling, dependencies);
end

%G %Display the content of matrix G

%% Fill the F matrix with capacitors from netlist

dep = namedargs2cell(dependencies);

% allowed_deps = {'xpos', 'vihc_ss', 'vohc_ss'};
allowed_deps = {'xpos'}; % elements depending on vihc_ss, vohc_ss need to be evaluated later

for i = 1:numCap
    
    if omit_state_dep_elements == true
        
        deps_ok = true;
        % allow only resistors that depend on nothing or xpos
        for j = 1:numel(Capacitor(i).dep)
            % if all(strcmp(Capacitor(i).dep(j), allowed_deps) == 0)
            %     deps_ok = false;
            % else
            if is_in(Capacitor(i).dep{j}, allowed_deps)
            else
                deps_ok = false;
            end
        end

        if deps_ok == false
            continue
        end
    end
    
    n1 = Capacitor(i).node1;
    n2 = Capacitor(i).node2;
    
    f = Capacitor(i).value(dep{:});
    
    f = f * NumstacksScaling; % parallel capacitors
    
    %If neither side of the element is connected to ground
    %then subtract it from appropriate location in matrix.
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

%% Fill the I matrix
% [I{:}]=deal('0');
for j = 1:numNode
    for i = 1:numI
        if Isource(i).node1 == j
            I(j) = I(j) - Isource(i).value(dep{:});
        elseif Isource(i).node2 == j
            I(j) = I(j) + Isource(i).value(dep{:});
        end
    end
end

%% Fill the V matrix
% for i = 1:numNode
%       V(i) = Vsource(i).Value;
% %     don't need this???
% end

%%
if numV ~= 0
    
    %% Fill the B matrix
    %First handle the case of the independent voltage sources.
    for i = 1:numV            %Go through each independent source.
        for j = 1:numNode     %Go through each node.
            if Vsource(i).node1 == j        %If node is first node,
                B(j,i) = 1;                 %then put '1' in the matrices.
            elseif Vsource(i).node2 == j    %If second node, put -1.
                B(j,i) = -1;
            end
        end
    end
    
    %% Fill the C matrix
    % %First handle the case of the independent voltage sources.
    % for i = 1:numV            %Go through each independent source.
    %     for j = 1:numNode     %Go through each node.
    %         if Vsource(i).node1 == j       %If node is first node,
    %             C(i,j) = 1;                %then put '1' in the matrices.
    %         elseif Vsource(i).node2 == j   %If second node, put -1.
    %             C(i,j) = -1;
    %         end
    %     end
    % end
    % In this setup actually a transpose of B
    C = B';
    
    %% Fill the D matrix
    %The D matrix is non-zero only for CCVS and VCVS (not included
    %in this simple implementation of SPICE)
    
    %% Fill the E matrix
    for i=1:numV
        E(i) = Vsource(i).value(dep{:});
    end
    
    %% Fill the H matrix
%     for i=1:numV,
%         H(i)=I(i); %Isource(i).Value;
%     end

end

%% Form the A, and Z matrices (As cell arrays of strings).

if numV ~= 0
    Acell = [G B; C D]; %a complete conductance matrix from resistors
    Gcapcell = [F J; K N]; %a complete conductance matrix from capacitors
    Zcell = [I; E];     %a complete matrix of independent sources
else
    Acell = G;
    Gcapcell = F;
    Zcell = I;
end

end
