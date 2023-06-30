function C = CapacitanceMatrix( t, y, Z, ...
        C0, channels, Numstacks, spec_size, ...
        c, Ic, xpos, dep_fun, p, q, mechopt)
% CAPACITANCEMATRIX Calculates MNA C matrix
% C dy/dt = -A(t,y)*y + z = f(x,t)

% inverse permutations
P(p) = 1:numel(p);
Q(q) = 1:numel(q);

% Extract IHC & OHC receptor potentials from solution
fn = fieldnames(dep_fun);
for i = 1:numel(fn)
    dep = fn{i};
    vars.(dep) = dep_fun.(dep)(y(Q));
end

% Link to steady-state solution
vars.vihc_ss = mechopt.vihc_ss;
vars.vohc_ss = mechopt.vohc_ss;

C = C0;

n = spec_size * Numstacks;

% deps = fieldnames(c);
% deps = {'vihc', 'vohc'};
deps = {};
possilble_deps = fieldnames(Ic);
for j = 1:numel(possilble_deps)
    dep = possilble_deps{j};
    if ~isempty(Ic.(dep))
        deps{end+1} = dep;
    end
end


if ~isempty(channels)
    deps = [deps, channels.var_name];

    % Fill in Voltage dependent part of C
    for j = 1:numel(deps)
        dep = deps{j};
        v = vars.(dep);
        assert(all(size(xpos) == size(v)));
        for i = 1:numel(c.(dep))
            cc = c.(dep){i}(v);
            assert(all(size(cc) == size(v)));
            C(1:n,1:n) = C(1:n,1:n) + kron(spdiags(cc,0,Numstacks,Numstacks),Ic.(dep){i});
        end
    end

    V = struct('IHC', vars.vihc, 'OHC', vars.vohc, 'IHC_ss', vars.vihc_ss, 'OHC_ss', vars.vohc_ss);

    % M is not 'Capacitance' anymore
    M = ChannelspOpen(V, channels, Numstacks, ...
            'mass', true, 'system', false, 'rhs', false, 'jacobian', false);

    [m,~] = size(M);

    C((n+1):(n+m),(n+1):(n+m)) = M;

    % C = sparse(blkdiag(C, mechopt.mass));
end

end