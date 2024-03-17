function N_Na=fxnint(markovrates,Na_max,dt)
%FXNINT 
%   Calculate number of open sodium channels according to Eqns. (9)-(13)
%   of Mino et al. (2002) with nint(N_Na)

global m h;

n = size(markovrates.am,1);

m_std = sqrt((2/Na_max) * (markovrates.am .* markovrates.bm) ./ (markovrates.am + markovrates.bm));
h_std = sqrt((2/Na_max) * (markovrates.ah .* markovrates.bh) ./ (markovrates.ah + markovrates.bh));

% for i = 1:n
%     mnew = -1;
%     while (mnew < 0) || (mnew > 1) %!!!
%         mnew = m(i) + ((markovrates.am(i) .* (1-m(i))) - (markovrates.bm(i) .* m(i)))*dt + m_std(i)*randn(1)*sqrt(dt);
%     end
%     m(i) = mnew;
% 
%     hnew = -1;
%     while (hnew < 0) || (hnew > 1) %!!!
%         hnew = h(i) + ((markovrates.ah(i) .* (1-h(i))) - (markovrates.bh(i) .* h(i)))*dt + h_std(i)*randn(1)*sqrt(dt);
%     end
%     h(i) = hnew;
% end

% VECTORIZED
mnew = -1 * ones(n,1);
ind = 1:n;
while ~isempty(ind) %any(mnew < 0) || any(mnew > 1) %!!!
    q = numel(ind);
    mnew(ind) = m(ind) + ((markovrates.am(ind) .* (1-m(ind))) - (markovrates.bm(ind) .* m(ind)))*dt + m_std(ind) .* randn(q,1) * sqrt(dt);
    
    ind = find( mnew < 0 | mnew > 1 );    
end
m = mnew;

hnew = -1 * ones(n,1);
ind = 1:n;
while ~isempty(ind) %any(hnew < 0) || any(hnew > 1) %!!!
    q = numel(ind);
    hnew(ind) = h(ind) + ((markovrates.ah(ind) .* (1-h(ind))) - (markovrates.bh(ind) .* h(ind)))*dt + h_std(ind) .* randn(q,1) * sqrt(dt);
    
    ind = find( hnew < 0 | hnew > 1 );
end
h = hnew;

% ORIGINAL !!!
% mnew = -1;
% while (mnew < 0) || (mnew > 1)
%     mnew = m + ((markovrates.am .* (1-m)) - (markovrates.bm .* m))*dt + m_std*randn(1)*sqrt(dt);
% end
% m = mnew;
% 
% hnew = -1;
% while (hnew < 0) || (hnew > 1)
%     hnew = h + ((markovrates.ah .* (1-h)) - (markovrates.bh .* h))*dt + h_std*randn(1)*sqrt(dt);
% end
% h = hnew;

N_Na = round(Na_max * m.^3 .* h);
find(N_Na > 0)

end

