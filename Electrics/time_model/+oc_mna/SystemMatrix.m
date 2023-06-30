function B = SystemMatrix(B0, Numstacks, G, I, xpos, vars, size_info)
% SYSTEMMATRIX Calculates MNA A matrix (mostly conductivity matrix)
% C dy/dt = -A(t,y)*y + z = f(x,t)
%         B := -A
%        B0 := -A0
arguments
    B0 (:,:) double
    Numstacks (1,1) double
    G (1,1) struct
    I (1,1) struct
    xpos (:,1) double
    vars (1,1) struct
    size_info (1,1) struct
end
% persistent spdiag_N
% if isempty(spdiag_N)
%     spdiag_N = spdiags(ones(Numstacks,1),0,Numstacks,Numstacks);
% end

B = B0;

deps = fieldnames(G);

m = size_info.circuit.start;
n = size_info.circuit.end;

% nnzB0 = nnz(B0);
nzmaxB0 = nzmax(B0);

deltaB = [];

for j = 1:numel(deps)
    dep = deps{j};
    v = vars.(dep);
    assert(all(size(xpos) == size(v)));
    % Fill in 'dep'-dependent part of B
    for i = 1:numel(G.(dep))
        g = G.(dep){i}(v);
        assert(all(size(g) == size(v)));
        
        % deltaA = kron(spdiag_N .* g, I.(dep){i}); % slower
        % deltaA = kron(spdiags(g,0,Numstacks,Numstacks),I.(dep){i});
        % B(m:n,m:n) = B(m:n,m:n) - deltaA;

        if isempty(deltaB)
            deltaB = kron(spdiags(g,0,Numstacks,Numstacks),I.(dep){i});
        else
            deltaB = deltaB + kron(spdiags(g,0,Numstacks,Numstacks),I.(dep){i});
        end
    end
end
B(m:n,m:n) = B(m:n,m:n) - deltaB;

assert(nzmax(B) == nzmaxB0)
end