function [G0S, G0M, hfig] = g0(x, plotflag)
%
% Computes coarse profiles of the hydrodynamic Green's function for a
% realistic geometry of the human cochlear duct.
%
% [G0S, G0M] = G0(X, PlotFlag)
%
% x   = BM point vector (points may be irregularly spaced!)
% plotflag = if set to 1, plots are shown [def.0]
%
% G0S = Stapes_BM  coarse Green's function (Stapes-BM  fluid coupling);
% G0M = Basilar membrane coarse Green's function (BM-BM fluid coupling).
%
% G0S and G0M are the mere result of the integration over the inverse of the
% coclear duct section surfaces and are NOT multiplied by the BM width.
% This latter step is accomplished in the routine GREENF.M. The dimension of
% G0S and G0M is  1/BETA  in BETA units (1 BETA = 33.5 mm).
%
% (R.Nobili-Padova University, F.Mammano-SISSA, rev. 21-11-97)

if nargin < 2
	plotflag = 0;
end

if nargin < 1
	plotflag = 1;
    x = (1:400)/400;
end

%% UNIT CONVERSION

BETA = 33.5; % [mm]

%%
% disp('	Computing G0 ...');

x = x(:)'; 	                % Make sure x is a row
N = length(x);

% yp, ym is cochlear duct radius at the positions xp, xm
[xp, yp, xm, ym, ap, am] = mech.cochlearduct();

% Find cross section areas at middle points from human cochlea data
% csa_switch = 2;
csa_switch = 1; % used by Vetesnik

if csa_switch == 1
    % assuming cochlear duct cross sections are circles in this version
    ypm = (yp(1:end-1) + yp(2:end))/2;
    Sp = pi*ypm.^2;

    ymm = (ym(1:end-1) + ym(2:end))/2;
    Sm = pi*ymm.^2;
else
    Sp = ap(2:end) - diff(ap)/2;
    Sm = am(2:end) - diff(am)/2;
end

dx_over_Sp = diff(xp)./Sp;		% Find integrand
dx_over_Sm = diff(xm)./Sm;

Intp_end = (0.5*pi)*yp(end)/Sp(end);	% Contribution from upper side of helicotrema
Intm_end = (0.5*pi)*ym(end)/Sm(end); 	% Contribution from lower side of helicotrema

% emulate a mistake (?)
% Intp_end = Intm_end;
% Intm_end = 0;

% integrate
Intp = cumsum([dx_over_Sp; Intp_end], 'reverse')';
Intm = cumsum([dx_over_Sm; Intm_end], 'reverse')';

% Interpolation of Intp and Intm:
%G0S = interpol(xp,Intp,x) + interpol(xm, Intm,x);

Fp = griddedInterpolant(xp, Intp);
Fm = griddedInterpolant(xm, Intm);

G0S = Fp(x) + Fm(x);

G0M = G0S'*ones(size(G0S));
G0M = tril(G0M) + tril(G0M)'- diag(G0S);

if plotflag == 1
    hfig(1) = figure;
    hold on
    plot(Sp)
    plot(Sm)
    % plot(2*ap(2:end)-diff(ap), 'r:')
    % plot(2*am(2:end)-diff(am), 'b:')
    % plot(ap(2:end)-diff(ap)/2, 'r--')
    % plot(am(2:end)-diff(am)/2, 'b--')
    xlabel('Fractional distance from stapes');
	ylabel('');
	title('Cochlear duct crossection areas');
    
    hfig(2) = figure;
	plot(xp, yp*BETA, xm, -ym*BETA)
	xlabel('Fractional distance from stapes');
	ylabel('mm');
	title('Radii of human cochlea scalae');
    
    J = unique([1:20:N,N]);
    cmap = rwb(numel(J));
    
	hfig(3) = figure;
    hold on
    for j = 1:numel(J)
        plot(x', G0M(:,J(j)), 'Color', cmap(j,:), 'LineStyle', ':')
    end
	% plot(x', G0M(:,(1:20:N)),'b:')
	title('Coarse BM-BM pressure Green''s functions'),
    
    [b, ~, hfig(4)] = mech.bmwidth(x, plotflag);
    b = b/2;  % Effective BM width in BETA units

    Gc = G0M.*(b'*b);
	hfig(5) = figure;
    hold on
    for j = 1:numel(J)
        plot(x', Gc(:,J(j)), 'Color', cmap(j,:))
    end
	% plot(x',Gc(:,1:20:N)),
    title('Coarse BM-BM force Green''s function')
else
    hfig = [];
end

end