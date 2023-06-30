function [ I ] = Gt_pattern( Resistor, numNode, nn, mm )
%GT_PATTERN
% I containts a pattern where time-dependent values of circuit are stored
% in the matrix. It contains a sign, but not value. The third dimmension of
% I corresponds to the number of time dependent elements.


% extract time-dependent resistors
% Resistor_t_dep = Resistor([Resistor.time_dep] == true);
Resistor_t_dep = Resistor;

n = numel(Resistor_t_dep);

% I = zeros(numNode,numNode,n);
I = zeros(nn,mm,n);

for i = 1:n
%     if Resistor_t_dep(i).time_dep == true
        
        n1 = Resistor_t_dep(i).node1;
        n2 = Resistor_t_dep(i).node2;
        
        % If neither side of the element is connected to ground
        % then subtract it from appropriate location in matrix.
        if (n1 ~= 0) && (n2 ~= 0)
            I(n1,n2,i) = -1;
            I(n2,n1,i) = -1;
        end

        % If node 1 is connected to ground, add element to diagonal
        % of matrix.
        if (n1 ~= 0)
            I(n1,n1,i) = 1;
        end
        % Ditto for node 2.
        if (n2 ~= 0)
            I(n2,n2,i) = 1;
        end
%     end
end

end
