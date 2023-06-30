function jac = Jacob(t, y, ...
        z, J0, channels, Numstacks, size_info, mechfs, ...
        G, dG, I, xpos, dep_fun, dep_sel, dep_sign, ...
        dep_names, dep_names_alt, ...
        p, q, needs_ordering, ...
        mechopt, t_mech, loaded_BM_displ, loaded_OHC_cilia)
arguments
    t (1,1) double
    y (:,1) double
    z (:,1) double
    J0 (:,:) double
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
% JACOB Calculates Jacobian
%
%  J = df/dy =
%    [ df1/dy1  df1/dy2  df1/dy3  ... df1/dyn ]
%    [ df2/dy1  df2/dy2  df2/dy3  ... df2/dyn ]
%    [ ...                                    ]
%    [ dfm/dy1  dfm/dy2  dfm/dy3  ... df2/dyn ]
%
% where
%  C dy/dt = -A(t,y)*y + z(t,y) = f(x,t) .
%
% If A(t,y) = A(t)
%  J = -A(t) + dz(t,y)/dy
%
% Otherwise: product rule
%  J = -A(t,y) - dA(t,y)/dy * y + dz(t,y)/dy
%
% Note: h(y) = -A*y is matrix * vector, i.e.
%
% h_i = - A_{i,1} y_1 - A_{i,2} y_2 - ... - A_{i,n} y_n 
%     = - sum_j^n A_{i,j} y_j
%
%    (dh/dy)_{i,k} = dh_i/dy_k = - A_{i,k} - sum_j^n dA_{i,j}/dy_k y_j
%                              = - A_{i,k} - dA_{i,:}/dy_k * y
%
% but often we know that for most j: dA_{i,j}/dy_k = 0
%
% so (dh/dy)_{i,k} = dh_i/dy_k = - A_{i,k} - sum_{j \in K} dA_{i,j}/dy_k y_j
%

% spec_size = mnaopt.circuit.num_v + mnaopt.circuit.num_node
spec_size = size_info.circuit.spec_size;

%%

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
        [BM_displ, OHC_cilia] = BMdispl( ...
            t, t_mech, mechfs, loaded_BM_displ, loaded_OHC_cilia);

    case 'electric'
        [BM_displ, OHC_cilia] = mech_interpolate( ...
            t, y, xpos, size_info, mechopt);
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

vars.xpos = xpos;

args_names = fieldnames(vars);
args = struct2cell(vars);

%% MNA equations: A(t,y)

% J0 contains constant part of jacobian and the minus (!) -A0 matrix
% i.e. J0 = -A0 + Jconst
% The function SystemMatrix subtracts the time and y-dependent parts of A
% i.e. jac = J0 - Delta A(t,y)

jac = oc_mna.SystemMatrix(J0, Numstacks, G, I, xpos, vars, size_info);

nzmax_orig = nzmax(jac);

% This is an older alterinative version
% A = oc_mna.SystemMatrix(1, A0, Numstacks, G, I, xpos, vars, size_info);
% jac = spalloc(size(A,1), size(A,2), nzmax(A));
% jac(1:end, 1:end) = -A;
% assert(nzmax(jac) == nzmax(A))

%% MNA equations: dA/dy * y

% V4
% let g(x;H,S) = boltzman(x,1,H,S);
% and x = v - w
% then
% dg/dv = boltzman_der(v,1,H+w,S) = boltzman_der(x,1,H,S)
% and
% dg/dw = -boltzman_der(-w,1,H-v,S) = -boltzman_der(x,1,H,S)
%
% therefore we can modify V1:
%

mm = size_info.circuit.start;
nn = size_info.circuit.end;
kk = size_info.circuit.size;

deps = {'vihc', 'vohc'};

for k = 1:numel(deps)
    dep = deps{k};

    n_elem = numel(I.(dep));
    if n_elem == 0
        continue
    end
    
    assert(nzmax(jac) - nnz(jac) >= 4*Numstacks); % 4 should be something like numel(TT) * n_elem
    
    s = dep_sign.(dep);
    w = [dep_sel.(dep){1}(1), dep_sel.(dep){2}(1)];
    
    for i = 1:n_elem
        
        TT = I.(dep){i};
        dGG = dG.(dep){i};

        v = dGG.nodes_affected;

        [~,j] = cellfun(@(x) ismember(x, args_names), dGG.args , 'UniformOutput', true);

        dg = dGG.fun(args{j});
        
        J = spalloc(kk,kk,Numstacks*numel(v)*2);

        for n = 1:Numstacks

            l = (n-1)*spec_size; % ~ which cross-section

            for vv = v
                tmp = dg(n) * TT(vv,v) * y(v+l);
                J(l+vv,l+w(1)) = s(1) * tmp; % s(1) was previously 1
                J(l+vv,l+w(2)) = s(2) * tmp; % s(2) was previously -1
            end
        end    

        jac(mm:nn,mm:nn) = jac(mm:nn,mm:nn) - J;
    end
end

%% Dependence of circuit (MNA equations) on
%    1. The IHC/OHC basolateral open probability (channels block)
%    2. The BM displacement and OHC stereocilia (mech block)

mm = size_info.circuit.start;
nn = size_info.circuit.end;
kk = size_info.circuit.size;

deps_consts = struct();

deps = dep_names;
deps_alt_names = dep_names_alt;

deps_consts.popen_IHC_fast = 1;
deps_consts.popen_IHC_slow = 1;
deps_consts.popen_IHC_n = 1;
deps_consts.popen_IHC_slow_1 = 1;
deps_consts.popen_IHC_slow_2 = 1;
% deps_consts.popen_OHC_fast = 1;
% deps_consts.popen_OHC_slow = 1;
deps_consts.popen_OHC_A = 1;
deps_consts.popen_OHC_C = 1;

if strcmp(mechopt.integration, 'electric') % && ~any(isnan(mechopt.vohc_0))
    
    deps = [deps, ...
        {'BM_displ', 'OHC_cilia'}];

    deps_alt_names = [deps_alt_names, ...
        {'mech_BMx', 'mech_TMx'}];

    deps_consts.BM_displ = mechopt.NMkonst_BM;
    deps_consts.OHC_cilia = mechopt.NMkonst_OHC_cilia;

end

for k = 1:numel(deps)
    dep = deps{k};
    
    n_elem = numel(I.(dep));
    if n_elem == 0
        continue
    end

    mmm = size_info.(deps_alt_names{k}).start;
    nnn = size_info.(deps_alt_names{k}).end;
    kkk = size_info.(deps_alt_names{k}).size;

    assert(nzmax(jac) - nnz(jac) >= 4*kkk); % 4 should be something like numel(TT) * n_elem

    K = deps_consts.(dep); % this is basically d\xi/dy
    
    sel = dep_sel.(dep);
    sel = sel - mmm + 1;

    for i = 1:n_elem

        TT = I.(dep){i};
        dGG = dG.(dep){i};

        [~,j] = cellfun(@(x) ismember(x, args_names), dGG.args , 'UniformOutput', true);

        dg = K * dGG.fun(args{j});
        v = dGG.nodes_affected;

        ver = 3; % version selector
        if ver == 1
            % original version (slow)
            J = spalloc(kk,kkk,Numstacks*numel(v)*2);        
            for n = 1:Numstacks
                l = (n-1)*spec_size; % ~ which cross-section
                for vv = v
                    jj = dg(n) * TT(vv,v) * y(v+l);
                    J(l+vv,sel(n)) = jj;
                end
                % the above actually expands to ... which is even slower
                % for vv = v
                %     jj = 0;
                %     for vvv = v
                %         jj = jj + dg(n) * TT(vv,vvv) * y(vvv+l);
                %     end
                %     J(l+vv,sel(n)) = jj;
                % end
            end
            
        elseif ver == 2
            % improved version (slightly faster)
            V = zeros(1,Numstacks*numel(v));
            i1 = zeros(1,Numstacks*numel(v));
            i2 = [sel; sel];
    
            for n = 1:Numstacks
                l = (n-1)*spec_size; % ~ which cross-section
                for o = 1:numel(v)
                    vv = v(o);
                    jj = dg(n) * TT(vv,v) * y(v+l);
                    V((n-1)*numel(v) + o) = jj;
                    i1((n-1)*numel(v) + o) = l + vv;
                end
            end
            J = sparse(i1, i2(:), V(:), kk, kkk);

        elseif ver == 3
            % vectorized version (fast)
            assert(iscolumn(dg))
            assert(iscolumn(y))
            assert(isrow(sel))

            i1 = v(:) + ((1:Numstacks)-1)*spec_size;
            i2 = [sel; sel];
    
            VV = zeros(numel(v), Numstacks);
            for o = 1:numel(v)
                vv = v(o);
                VV(o,:) = dg' .* (TT(vv,v) * y(i1));
            end
            VV = VV(:);
            J = sparse(i1(:), i2(:), VV, kk, kkk);
        end

        jac(mm:nn,mmm:nnn) = jac(mm:nn,mmm:nnn) - J;
        

    end
end

%%

% BACKUP - this is now included in the block above
% switch mechopt.integration
%     case 'standalone'
%     case 'electric'
%         if ~isnan(mechopt.vohc_0)
%             
%             % spec_sel = (0:(Numstacks-1)) * spec_size;
% 
%             mm = size_info.circuit.start;
%             nn = size_info.circuit.end;
%             kk = size_info.circuit.size;
% 
%             deps = {'BM_displ', 'OHC_cilia'};
%             deps_alt_names = {'mech_BMx', 'mech_TMx'};
%             
%             deps_consts = struct( ...
%                 'BM_displ', mechopt.NMkonst_BM, ...
%                 'OHC_cilia', mechopt.NMkonst_OHC_cilia);
% 
%             for k = 1:numel(deps)
%                 dep = deps{k};
%                 
%                 n_elem = numel(I.(dep));
%                 if n_elem == 0
%                     continue
%                 end
% 
%                 mmm = size_info.(deps_alt_names{k}).start;
%                 nnn = size_info.(deps_alt_names{k}).end;
%                 kkk = size_info.(deps_alt_names{k}).size;
% 
%                 assert(nzmax(jac) - nnz(jac) >= 4*kkk); % 4 should be something like numel(TT) * n_elem
% 
%                 K = deps_consts.(dep); % this is basically d\xi/dy
%                 
%                 sel = dep_sel.(dep);
%                 sel = sel - mmm + 1;
% 
%                 for i = 1:n_elem
% 
%                     TT = I.(dep){i};
%                     dGG = dG.(dep){i};
% 
%                     v = dGG.nodes_affected;
% 
%                     [~,j] = cellfun(@(x) ismember(x, args_names), dGG.args , 'UniformOutput', true);
% 
%                     dg = K * dGG.fun(args{j});
% 
%                     J = spalloc(kk,kkk,Numstacks*numel(v)*2);
% 
%                     for n = 1:Numstacks
% 
%                         l = (n-1)*spec_size; % ~ which cross-section
% 
%                         for vv = v
%                             jj = dg(n) * TT(vv,v) * y(v+l);
%                             J(l+vv,sel(n)) = jj;
%                         end
%     %                     for vv = v
%     %                         jj = 0;
%     %                         for vvv = v
%     %                             jj = jj + dg(n) * TT(vv,vvv) * y(vvv+l);
%     %                         end
%     %                         jac(l+vv,sel(n)) = jj;
%     %                     end
%                     end
% 
%                     jac(mm:nn,mmm:nnn) = jac(mm:nn,mmm:nnn) - J;
% 
%     %                 J = spalloc(kk,kkk,Numstacks*numel(v)*2);
%     %                 
%     %                 for vv = v
%     %                     jj = zeros(Numstacks,1);
%     %                     for vvv = v
%     %                         jj = jj + dg .* TT(vv,vvv) .* y(vvv+spec_sel);
%     %                     end
%     %                     ind = sub2ind(spec_sel+vv,sel);
%     %                     J(ind) = jj(:)';
%     %                 end
%     %                 
%     %                 jac(mm:nn,mmm:nnn) = jac(mm:nn,mmm:nnn) - J;
% 
%                 end
%             end
%         end
%     otherwise
%         error('Unsupported value of mechopt.integration = %s', mechopt.integration)
% end

%% Extension
if ~isempty(channels)
    orders = [channels.order];

    V = struct('IHC', vars.vihc, 'OHC', vars.vohc, 'IHC_ss', vars.vihc_ss, 'OHC_ss', vars.vohc_ss);

    for i = 1:numel(channels)

        [~, ~, ~, JJ, dP_V, dZ_V, ~] = pOpenIHC(V.(channels(i).voltage), channels(i), ...
            'mass', false, 'system', false, 'rhs', false, 'jacobian', true);

        a = (size_info.channels.start - 1) + sum(orders(1:i-1)) * Numstacks;
        b = (size_info.channels.start - 1) + sum(orders(1:i-1)) * Numstacks;

        [aa,bb] = size(JJ);

        % dPdVfull = sparse(aa,b);

        m = Numstacks * sum(orders(1:i-1));
        ll = (size_info.channels.start - 1) + m;


        assert(channels(i).order == 2)

        popen_is_constant = false;
        switch channels(i).voltage
            case {'IHC_ss', 'OHC_ss'}
                popen_is_constant = true;
            case 'IHC'
                s = dep_sign.vihc;
                w = [dep_sel.vihc{1}(1), dep_sel.vihc{2}(1)];
            case 'OHC'
                s = dep_sign.vohc;
                w = [dep_sel.vohc{1}(1), dep_sel.vohc{2}(1)];
            otherwise
                error('Not set up!')
        end

        if ~popen_is_constant

            %v = dGG.nodes_affected;
            v = [2]; % only second eq. affected every cross-section
    
            dPdV = spalloc(aa,b,Numstacks*numel(v)*2);
    
            for n = 1:Numstacks
    
                l = (n-1)*spec_size; % ~ which cross-section
    
                for vv = v % only second eq. affected every cross-section
                    tmp = (dZ_V(n) + dP_V(n) * y(ll+n*vv));
                    dPdV(n*vv,l+w(1)) = s(1) * tmp;
                    dPdV(n*vv,l+w(2)) = s(2) * tmp;
                end
            end    
    
            % dPdVfull = dPdVfull + dPdV;
            jac((a+1):(a+aa),1:b) = jac((a+1):(a+aa),1:b) + dPdV;
            
        end

        jac((a+1):(a+aa),(b+1):(b+bb)) = jac((a+1):(a+aa),(b+1):(b+bb)) + JJ;    
        % jac =  [      jac, sparse(a,bb);
        %          dPdVfull,          JJ];

    end
end
%% Mechanical model jacobian
% OHC aplifier - depends on VOHC
%
% constant part of mech jacobian already in A0 and J0

if strcmp(mechopt.integration, 'electric') ...
    && strcmp(mechopt.approximation, 'nonlinear')

    switch mechopt.amplifier
        case 'none'
        case 'mechanic'            
            mm = size_info.mech.start;
            nn = size_info.mech.end;
            
            Jnl = mechopt.jacobian_nonlin_part(y(mm:nn));

            jac(mm:nn,mm:nn) = jac(mm:nn,mm:nn) + Jnl;
        case 'electric'

            if isnan(mechopt.vohc_0)
                eohc = zeros(size(vars.vohc));
            else
                eohc = vars.vohc - mechopt.vohc_0;
            end

            if ~isnan(mechopt.vohc_0)
            
                [J2, J4] = mechopt.jacobian_nonlin_part_vohc(eohc);

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
                                + spdiags(s(k) * J(:), 0, kkk, kkk);% .* dep_sel.vohc{k}(:);
                    end
            
                end
            end
        otherwise
            error('Unexpected value')
    end
end

%% Extend system by regularization term

if true%false
    
    w = [dep_sel.vohc{1}(1), dep_sel.vohc{2}(1)];
    s = dep_sign.vohc;

    jac(1,w(1)) = s(1) * (vars.vohc(1) > -50e-3);
    jac(1,w(2)) = s(2) * (vars.vohc(1) > -50e-3);

    s = size_info.reg.start:size_info.reg.end;
    jac(s,s) = 1;
end

%% Ordering

if needs_ordering
    jac = jac(p,q);
end


if any(isnan(jac(:)))
  error('nan encountered')
end

assert(nzmax_orig == nzmax(jac))

%%
plot_jac_spy = false;
if plot_jac_spy
    
    figure;
    hold on

    blocks = size_info.blocks;

    cmap = parula(numel(blocks)^2);
    cmap = jet(numel(blocks)^2);

    color = {
        [1,0,0],[1,1,0],[0,0.3,1],
        }

    k = 0;
    for i = 1:numel(blocks)
        b = blocks{i};
        m1 = size_info.(b).start;
        m2 = size_info.(b).end;

        for j = 1:numel(blocks)
            c = blocks{j};
            n1 = size_info.(c).start;
            n2 = size_info.(c).end;

            x = [m1 m2 m2 m1]
            y = [n1 n1 n2 n2]

            k = k + 1;
            cc = (color{i} + color{j})/2;
            a = 0.3
            if i >= j
                cc = (color{i} + color{j})/2;
            else
                cc = ([1,1,1] - color{i} + [1,1,1] - color{j})/3;
            end
            patch(x, y, cc, ...
                'EdgeColor', 'none', ...
                'FaceAlpha', a, ...
                'DisplayName', sprintf('%d: d_{%s} / d_{%s}', k, c, b))

        end
    end    

    spy(jac, 'k.', 3)

    legend('Location','eastoutside')
    
    ax = gca()
    ax.XTick = []
    ax.YTick = []
    ax.XAxis.Visible='off'
    ax.YAxis.Visible='off'
end

%% V3
% let g(x;H,S) = boltzman(x,1,H,S);
% and x = v + w
% then
% dg/gv = boltzman_der(v,1,H-w,S) = boltzman_der(x,1,H,S)
% as well as
% dg/gw = boltzman_der(w,1,H-v,S) = boltzman_der(x,1,H,S)
%
% therefore we can modify V1:
%

% args = {xpos, VoltageIHC, VoltageOHC};
% args_names = {'xpos', 'vihc', 'vohc'};
% 
% w_ihc = [sel_ihc{1}(1), sel_ihc{2}(1)];
% w_ohc = [sel_ohc{1}(1), sel_ohc{2}(1)];
% 
% T = [I_vihc, I_vohc];
% dG = [dg_vihc, dg_vohc];
% 
% for k = 1:numel(T)
% 
%     v = zeros(1,size(T{k},1));
%     for i = 1:size(T{k},1)
%         v(i) = any(T{k}(i,:));
%     end   
%     v = find(v);
% 
%     m = numel(y)/Numstacks;
% 
%     if any(strcmp(dG{k}.args, 'vihc'))
%         w = w_ihc;
%     elseif any(strcmp(dG{k}.args, 'vohc'))
%         w = w_ohc;
%     end
%     
%     [~,j] = cellfun(@(x) ismember(x, args_names), dG{k}.args , 'UniformOutput', true);
% 
%     J = spalloc(size(jac,1),size(jac,2),Numstacks*numel(v)*2);
%     
%     dg = dG{k}.fun(args{j});
%     
%     for n = 1:Numstacks
%         J((n-1)*m+v,(n-1)*m+w) = dg(n) * T{k}(v,v) * repmat(y(v+(n-1)*m),1,numel(w));
%     end
%     
%     jac = jac - J;
%     
% end
% 
% 
% jac = jac(p,q);

%% V2
% let g(x;H,S) = boltzman(x,1,H,S);
% and x = v + w
% then
% dg/gv = boltzman_der(v,1,H-w,S) = boltzman_der(x,1,H,S)
% as well as
% dg/gw = boltzman_der(w,1,H-v,S) = boltzman_der(x,1,H,S)
%
% therefore we can modify V1:
%

% if ~isempty(I_vihc)
%     T = I_vihc{1};
% 
%     v = zeros(1,size(T,1));
%     for i = 1:size(T,1)
%         v(i) = any(T(i,:));
%     end   
%     v = find(v);
% 
%     m = numel(y)/Numstacks;
% 
%     w = [sel_ihc_ic(1), sel_ihc_ec(1)];
% 
%     Ji = spalloc(size(jac,1),size(jac,2),Numstacks*numel(v)*2);
%     for n = 1:Numstacks
%         Ji((n-1)*m+v,(n-1)*m+w) = dg_vihc{1}(VoltageIHC(n)) * T(v,v) * repmat(y(v+(n-1)*m),1,numel(w));
%     %     Ji((n-1)*m+v,(n-1)*m+w) = boltzmann_der(Voltage(n),580e-9,-26,11) * T(v,v) * repmat(y(v+(n-1)*m),1,numel(w));
%     end
%     
%     jac = jac - Ji;
%     
% end
% 
% if ~isempty(I_vohc)
%     T = I_vohc{1};
% 
%     v = zeros(1,size(T,1));
%     for i = 1:size(T,1)
%         v(i) = any(T(i,:));
%     end   
%     v = find(v);
% 
%     m = numel(y)/Numstacks;
% 
%     w = [sel_ohc_ic(1), sel_ohc_ec(1)];
% 
%     Jo = spalloc(size(jac,1),size(jac,2),Numstacks*numel(v)*2);
%     for n = 1:Numstacks
%         Jo((n-1)*m+v,(n-1)*m+w) = dg_vohc{1}(xpos(n), VoltageOHC(n)) * T(v,v) * repmat(y(v+(n-1)*m),1,numel(w));
%     end
% 
%     jac = jac - Jo;
% end

%% V1

% JJ = J;
% 
% T = I{2};
% 
% v = zeros(1,size(T,1));
% for i = 1:size(T,1)
%     v(i) = any(T(i,:));
% end   
% v = find(v);
% 
% m = numel(y)/Numstacks;
% 
% J = spalloc(size(jac,1),size(jac,2),Numstacks*numel(v)*2);
% for n = 1:Numstacks
%     sel = (1+(n-1)*m) : n*m;
%     V_ihc_ic = y(sel_ihc_ic(n));
%     V_ihc_ec = y(sel_ihc_ec(n));
%     U = spalloc(size(T,1),size(T,2),numel(v)*2);
%     for vv = v
%         U(vv,sel_ihc_ic(1)) = (T(vv,v) .* boltzmann_der(V_ihc_ic,580e-9,-26-V_ihc_ec,11)) * y(v+(n-1)*m);
%         U(vv,sel_ihc_ec(1)) = (T(vv,v) .* boltzmann_der(V_ihc_ec,580e-9,-26-V_ihc_ic,11)) * y(v+(n-1)*m);
%     end
%     J(sel, sel) = U;
% end

%%

% J = JJ;
% 
% [m,n] = size(J);
% 
% J = zeros(size(jac));
% for j = 1:n % over columns
%     for i = 1:m % over rows
%         tmp = 0;
%         for r = 1:n
% %             tmp = tmp + dA(i,r)dy(j) * y(r)
%             tmp = tmp + dAirdyj(i,r,j) * y(r);
%         end
%         J(i,j) = tmp;
%     end
% end
% 
%     function d = dAirdyj(i,r,j)
%         if ~find(j == w)
%             d = 0;
%         elseif ~find(i == v)
%             d = 0;
%         elseif ~find(r == v)
%             d = 0;
%         else
%             d = T
%         end
%     end

%%
% this is not right... the partial jacobian should be block diagonal

% T = kron(spdiags(ones(size(gg)),0,Numstacks,Numstacks),I{2});
% 
% J = zeros(size(jac));
% for j = 1:size(J,2) % over columns
%     [~, k] = find(j == sel_ihc_ic, 1);
%     [~, l] = find(j == sel_ihc_ec, 1);
%     if k
%         V_ihc_ic = y(sel_ihc_ic(k));
%         V_ihc_ec = y(sel_ihc_ec(k));
%         for i = 1:size(J,1) % over rows
%             J(i,j) = (T(i,:) .* boltzmann_der(V_ihc_ic,580e-9,-26-V_ihc_ec,11)) * y;
%         end
%     end
%     if l
%         V_ihc_ic = y(sel_ihc_ic(l));
%         V_ihc_ec = y(sel_ihc_ec(l));
%         for i = 1:size(J,1) % over rows
%             J(i,j) = (T(i,:) .* boltzmann_der(V_ihc_ec,580e-9,-26-V_ihc_ic,11)) * y;
%         end
%     end
% end
% J = sparse(J);
%%

% jac = jac - J;

end
