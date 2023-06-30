function [MM, AA, zz, JJ, dPP_V, dZZ_V, MM_V] = ChannelspOpen(V, channels, Numstacks, do)
%CHANNELSPOPEN
arguments
    V struct
    channels struct
    Numstacks (1,1) double
    do.mass (1,1) logical = true
    do.system (1,1) logical = true
    do.rhs (1,1) logical = true
    do.jacobian (1,1) logical = true
end

nv = namedargs2cell(do);

% [ISTILDE,CALLERTXT,CALLER,COMPONENTS] = detectOutputSuppression(nargout)

orders = [channels.order];

total_order = sum(orders);

d = Numstacks * total_order;
dd = Numstacks * numel(channels);

[MM, AA, JJ, MM_V] = deal(sparse(d,d));
[zz] = deal(zeros(d,1));
[dPP_V, dZZ_V] = deal(zeros(dd,1));

for i = 1:numel(channels)

    
    [M, A, z, J, dP_V, dZ_V, M_V] = pOpenIHC(V.(channels(i).voltage), channels(i), nv{:});

    m = Numstacks * sum(orders(1:i-1)) + 1;
    n = m + Numstacks * orders(i) - 1;
    o = Numstacks * (i-1) + 1;
    p = o + Numstacks - 1;
    
    if do.system
        AA(m:n,m:n) = A;
    end
    if do.mass
        MM(m:n,m:n) = M;
        MM_V(m:n,m:n) = M_V;
    end
    if do.jacobian
        JJ(m:n,m:n) = J;
        dPP_V(o:p) = dP_V;
        dZZ_V(o:p) = dZ_V;
    end
    if do.rhs    
        zz(m:n) = z;
    end
    
    
end
end

