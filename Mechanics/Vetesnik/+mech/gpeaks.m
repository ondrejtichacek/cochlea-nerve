function [P, dx, hfig] = gpeaks(x, plotflag)
%
% 	PEAKS.M  computes the Green's function peak matrix of the human
%       cochlea. Data is loaded from MAN_DATA.MAT.
%
%		[P, dx] = GPEAKS(x);
%
% 	x  =  BM point vector.
%	plotflag = If set to 1, plots are showed (def. 0)
%
%	dx = Interpoint distances
%	P  = Peak matrix
%
%	Matrix P is not normalised (see PSIGMA.M) and must be multiplied
%	by b(x)*b(x'), where b(x) is the BM width. The dimension of P
%	is  1/BETA  in BETA units (1 BETA = 33.5 mm).
%
%	Files CDDATA.MAT, BMWDATA.MAT, INTERPOL.M  and PSIGMA.M are called.
%
%			(R.Nobili-Padova University, last rev. 15-2-03)

if nargin < 2
	plotflag = 0;
end

if nargin < 1
	plotflag = 1;
    x = (1:400)/400;
end

% disp('	Function GPEAKS is at work');

x = x(:);		% To make sure that x is a column
N = length(x);

c = [x(1);x(2:N)+x(1:N-1)]/2;	% averaged coordinates

b = mech.bmwidth(x);

% On account of radial BM curvature, the effective BM width is assumed geometric-width/2
b = b/2;

[xp, yp, xm, ym] = mech.cochlearduct();

% upper cochlear radius
% rp = interpol(xp, yp, c);
Fp = griddedInterpolant(xp, yp);
rp = Fp(c);

% lower cochlear radius
% rm = interpol(xm, ym, c);
Fm = griddedInterpolant(xm, ym);
rm = Fm(c);

r = (rp+rm)/2;   % Computing averaged coclear radius at point c(i)

xo = [0; x(1:N-1)];
dx = x-xo;	% interval lengths

[P, Q] = deal(zeros(N,N));

for j = 1:N
    cj = c(j);
    bj = b(j);
    rj = r(j);
    
    % The following routine computes peaks without reflection effect at the stapes
    % Q(:,j) = psigma(x-cj,bj,rj) - psigma(xo-cj,bj,rj);
    
    % The following routine computes peaks with reflection at stapes,
    % i.e.  peak + image-peak  mirrored a x = 0.
    P(:,j) = mech.psigma(x-cj,bj,rj) - mech.psigma(xo-cj,bj,rj) ...
           + mech.psigma(x+cj,bj,rj) - mech.psigma(xo+cj,bj,rj);
end

if plotflag == 1
    J = unique([1:20:N,N]);
	hfig = figure();
    hold on
	cmap = rwb(numel(J));    
    for j = 1:numel(J)
        plot(x, P(:,J(j)), 'Color', cmap(j,:))
    end
	title('BM-BM Green''s function peaks')
else
    hfig = [];
end

end