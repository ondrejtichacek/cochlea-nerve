function [ C, C0, G, G0, Z, Y0 ] = circuit_matrices( ...
        xgrid, Resistor, Capacitor, Vsource, Isource, ...
        numRes, numCap, numV, numI, numNode, ...
        Numstacks, numstacks_scaling, ...
        dependencies, sparse_flag )
%CIRCUIT_MATRICES calculates SPARSE steady-state matrix C of
% capacitors, G and G0 of resistors and SPARSE vectors Y0 of
% initial voltages V0 and Z of voltage sources in the circuit
%
% RETURNS:
% C : full capacitor matrix 
% G : full coupled resistor matrix
% G0: coupled resistor matrix without RAO and RAI (the varying elements)
% Z : full independent sources vector
% Y0: full initial voltages vector

arguments
    xgrid (1,:) double
    Resistor (1,:) CircuitElement
    Capacitor (1,:) CircuitElement
    Vsource (1,:) CircuitElement
    Isource (1,:) CircuitElement
    numRes (1,1) double
    numCap (1,1) double
    numV (1,1) double
    numI (1,1) double
    numNode (1,1) double
    Numstacks (1,1) double
    numstacks_scaling (1,1) double
    dependencies (:,1) struct
    sparse_flag (1,1) logical = true
end

% Create a cell for every cochlear cross section
[A, A0, C, C0, Z] = deal(cell(1,Numstacks));

for i = 1:Numstacks
    
    % a : conductance matrix from resistors
    % c : conductance matrix from capacitors
    % z : matrix of independent sources
    omit_state_dep_elements = false;
    [a, c, Z{i}] = RCconductances(Resistor, Capacitor, Vsource, Isource, ...
        numRes, numCap, numV, numI, numNode, omit_state_dep_elements, ...
        numstacks_scaling, dependencies(i));
    
    A{i} = sparse(a);
    C{i} = sparse(c);
    
    omit_state_dep_elements = true;
    [a0, c0, ~] = RCconductances(Resistor, Capacitor, Vsource, Isource, ...
        numRes, numCap, numV, numI, numNode, omit_state_dep_elements, ...
        numstacks_scaling, dependencies(i));
    
    A0{i} = sparse(a0);
    C0{i} = sparse(c0);
end

[n,~] = size(a);

% Matrices G, C are block-diagonal (without blocks' coupling). That means
% 
% G = [G{1} 0    0    0    ]
%     [0    G{2} 0    0    ]
%     [0    0    G{3} 0    ]
%     [0    0    0    G{4} ]
%
% where 0 is a zero matrix

G = blkdiag(A{:});
C = blkdiag(C{:});

G0 = blkdiag(A0{:});
C0 = blkdiag(C0{:});

Z = cat(1,Z{:});

dx = diff(xgrid);

%% now add links / coupling between blocks
% links are only formed for R_SV, R_SL and R_OC
if Numstacks > 1
    % num_long = sum(false == cellfun(@isempty, {Resistor.long_value}, 'UniformOutput', true));
    num_long = sum(true == cellfun(@(x) strcmp(x, 'longitudinal'), {Resistor.special}, 'UniformOutput', true));

    GG = spalloc(size(G,1), size(G,2), Numstacks*4*num_long);

    scaling = numstacks_scaling;

    for l = 1:numRes
        if strcmp(Resistor(l).special, 'longitudinal')
            k = Resistor(l).node1;
            % fprintf('Longitudinal connection of node %d\n', k);

            % the first Stack

            j = 1;

            scaling = dx(j) * Numstacks * numstacks_scaling;
            
            g = getResistorValue(Resistor(l), dependencies(j), scaling);
            
            p = k;
            GG(p,p) = GG(p,p) + g;

            for j = 2:(Numstacks-1)
                % add variable functions here *j; put in longitudinal values

                scaling = dx(j) * Numstacks;

                g = getResistorValue(Resistor(l), dependencies(j), scaling);
                
                p = n*(j-1) + k;
                q = n*(j-2) + k;

                GG(p,p) = GG(p,p) + g + g;
                GG(p,q) = GG(p,q) - g;
                GG(q,p) = GG(q,p) - g;
            end

            j = Numstacks-1;

            scaling = dx(j) * Numstacks;

            g = getResistorValue(Resistor(l), dependencies(j), scaling);
            
            % the last Stack        
            p = n*(j) + k;
            q = n*(j-1) + k;

            GG(p,p) = GG(p,p) + g;
            GG(p,q) = GG(p,q) - g;
            GG(q,p) = GG(q,p) - g;
        end
        
    end % of connectivity loop

    G0 = G0 + GG;

    G = G + GG;
end

%% Calculate Y0

Y0 = G\Z;

%%

if sparse_flag == false
    C = full(C);
    G = full(G);
    G0 = full(G0);
    Z = full(Z);
    Y0 = full(Y0);
end

%% ??????
% inconsistent initial conditions in node 1 in the 1th cross-section
% Use an inconsistent initial condition to test initialization.
% Spy0(1)=Spy0(1)- 0.0001*Spy0(1); %1e-4; 


    function g = getResistorValue(Res, dependencies, numstacks_scaling)
        
        dep = namedargs2cell(dependencies);
        
        R = Res.value(dep{:});

        switch Res.si_unit
            case 'S'
                g = R;
            case 'Ohm'
                g = 1/R;
            otherwise
                error('Unknown unit %s.', Res.si_unit);
        end

        g = g / numstacks_scaling; % serial resistors
    end

end