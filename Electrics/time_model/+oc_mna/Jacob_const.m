function jac = Jacob_const(Numstacks, size_info, dep_sel, dep_sign, p, q, mechopt)
arguments
    Numstacks (1,1) double
    size_info (1,1) struct
    dep_sel (1,1) struct
    dep_sign (1,1) struct
    p (1,:) double
    q (1,:) double
    mechopt (1,1) mechOpt
end
% JACOB Calculates constant part of Jacobian
%  J = df/dy =
%    [ df1/dy1  df1/dy2  df1/dy3  ... df1/dyn ]
%    [ df2/dy1  df2/dy2  df2/dy3  ... df2/dyn ]
%    [ ...                                    ]
%    [ dfm/dy1  dfm/dy2  dfm/dy3  ... df2/dyn ]
%
% where
%  C dy/dt = -A(t,y)*y + z(t,y) = f(x,t) .

% inverse permutations
P(p) = 1:numel(p);
Q(q) = 1:numel(q);

last_block = size_info.blocks(end);
N = size_info.(last_block).end;

jac = spalloc(N, N, 4*Numstacks);


%% Mechanical model jacobian

% OHC aplifier - depends on VOHC
if strcmp(mechopt.integration, 'electric') ...
    && strcmp(mechopt.amplifier, 'electric') ...
    && strcmp(mechopt.approximation, 'linear')
            
    [J2, J4] = mechopt.jacobian_lin_part_vohc();

    s = dep_sign.vohc;

    parts = {'mech_BMv', 'mech_TMv'};
    JJ = {J2, J4};

    assert(numel(s) == 2)
    for i = 1:numel(parts)
        part = parts{i};

        mmm = size_info.(part).start;
        nnn = size_info.(part).end;
        kkk = size_info.(part).size;
        
        J = JJ{i};

        for k = 1:numel(s)
            jac(mmm:nnn, dep_sel.vohc{k}) = jac(mmm:nnn, dep_sel.vohc{k}) ...
                    + spdiags(s(k) * J(:), 0, kkk, kkk);
        end

    end
end

if all(jac(:) == 0)
    jac = [];
end

end
