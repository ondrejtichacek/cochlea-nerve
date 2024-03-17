function [m, h, N_Na] = initcw(markovrates,Na_max,dt, hs)
%INITCW
%   Calculate initial number of open sodium channels according to Eqn. (17)
%   and Eqns. (A1)-(A6) of Mino et al. (2002)

% global Nmh;

n = size(markovrates.ah,1);

% Set initial number of channels in each state
if Na_max == 1e3 % Mean steady-state values known for Na_max = 1000
    Nmh = [246     6     0     0   731    17     0     0];
    settletime = 100; % Number of Euluer time steps over which to randomize initial conditions
else % otherwise determine approximately appropriate initial values
    % Calculate initial steady-state fraction of open h particles
    hinf = markovrates.ah ./ (markovrates.ah + markovrates.bh);        
    
    Nmh = zeros(n,8);
    Nmh(:,5) = round(Na_max*hinf);
    Nmh(:,1) = Na_max - Nmh(:,5);
    settletime = 200; % Number of Euluer time steps over which to randomize initial conditions
end

% Calculate transition rates of "escapes"
zeta = [3*markovrates.am + markovrates.ah,...
          markovrates.bm + markovrates.ah + 2*markovrates.am,...
        2*markovrates.bm + markovrates.ah + markovrates.am,...
        3*markovrates.bm + markovrates.ah,...
        3*markovrates.am + markovrates.bh,...
          markovrates.bm + markovrates.bh + 2*markovrates.am,...
        2*markovrates.bm + markovrates.bh + markovrates.am,...
        3*markovrates.bm + markovrates.bh];

% Vectors describing state transitions from current to next state
p_trans_curr = [0,1,2,2,3,3,4,1,5,2,6,3,7,4,8,5,6,6,7,7,8];
p_trans_next = [0,2,1,3,2,4,3,5,1,6,2,7,3,8,4,6,5,7,6,8,7];

t_l_total = zeros(1,n); % cumulative transition lifetime

Nmh_ensemble = Nmh;

% Randomize initial conditions
while (t_l_total <= settletime*dt)

    lambda = Nmh*zeta';

    t_l = - log(rand(1,n)) ./ lambda; % Determine transition lifetime

    % Calculate cumulative transition lifetime
    t_l_total = t_l_total + t_l; 

    % If cumulative transition lifetime is long enough, return initial
    % number of channels in each state
     if (t_l_total > settletime*dt)

        break;

    end

    % Calculate cumulative probability distribution for state
    % transitions
    P = cumsum([0,...
        3*markovrates.am .* Nmh(:,1),...
          markovrates.bm .* Nmh(:,2),...
        2*markovrates.am .* Nmh(:,2),...
        2*markovrates.bm .* Nmh(:,3),...
          markovrates.am .* Nmh(:,3),...
        3*markovrates.bm .* Nmh(:,4),...
          markovrates.ah .* Nmh(:,1),...
          markovrates.bh .* Nmh(:,5),...
          markovrates.ah .* Nmh(:,2),...
          markovrates.bh .* Nmh(:,6),...
          markovrates.ah .* Nmh(:,3),...
          markovrates.bh .* Nmh(:,7),...
          markovrates.ah .* Nmh(:,4),...
          markovrates.bh .* Nmh(:,8),...
        3*markovrates.am .* Nmh(:,5),...
          markovrates.bm .* Nmh(:,6),...
        2*markovrates.am .* Nmh(:,6),...
        2*markovrates.bm .* Nmh(:,7),...
          markovrates.am .* Nmh(:,7),...
        3*markovrates.bm .* Nmh(:,8)]/lambda);

    % Determine which state transition has occurred
    ind = find(rand(1) < P,1);

    Nmh(:,p_trans_curr(ind)) = Nmh(:,p_trans_curr(ind)) - 1; % decrease number of channels in current state
    Nmh(:,p_trans_next(ind)) = Nmh(:,p_trans_next(ind)) + 1; % increase number of channels in next state                     

    Nmh_ensemble = [Nmh_ensemble; Nmh]; % !!! NEED TO EDIT THIS !!!

end

N_Na = Nmh(:,end);

hs.Nmh = Nmh;

m = 0;
h = 0;


end

