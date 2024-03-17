function N_Na=initfxfloor(markovrates,Na_max,dt)
%INITFXFLOOR
%   Calculate initial number of open sodium channels according to Eqns.
%   (9)-(13) of Mino et al. (2002) with floor(N_Na)

    global m h;

    % Set initial values according to deterministic steady-state
    m = markovrates.am/(markovrates.am+markovrates.bm);
    h = markovrates.ah/(markovrates.ah+markovrates.bh);

    m_std = sqrt((2/Na_max)*(markovrates.am*markovrates.bm)/(markovrates.am+markovrates.bm));
    h_std = sqrt((2/Na_max)*(markovrates.ah*markovrates.bh)/(markovrates.ah+markovrates.bh));
    
    lps = 100;
       
    % Randomize initial conditions
    for lp=1:lps

        mnew = -1;
        while (mnew<0)||(mnew>1)
            mnew = m + ((markovrates.am*(1-m))-(markovrates.bm*m))*dt + m_std*randn(1)*sqrt(dt);
        end
        m = mnew;
        
        hnew = -1;
        while (hnew<0)||(hnew>1)
            hnew = h + ((markovrates.ah*(1-h))-(markovrates.bh*h))*dt + h_std*randn(1)*sqrt(dt);
        end
        h = hnew;
        
    end

    N_Na = floor(Na_max*m^3*h);
    
end

