function s = psigma(x, b, r)
%
%		S = PSIGMA(X, B, R)  (this is an adimensional function)
%
%	X = BM-location vector (function value is undetermined at X = 0);
%	B = BM-width vector;
%	R = cochlear radius vector.
%
%	S = Green's function peaks integrated over X and centered at X=0.
%
%	IF X1 and X2 are points on the BM, PSIGMA(X2,B,R) - PSIGMA(X1,B,R)
%	yields the area between X1 and X2 under the peak profile (centered
%	at X = 0).
%	Rectangular matrices X, B, R with same dimensions are accepted as
%	arguments of PSIGMA.
%	NOTE: the resulting matrix S should be multiplied by the form factor
%	matrix B(x)'*B(x), with B(x) the BM width vector ( ' = transposed).
%
%			(R.Nobili-Padova University, rev.21-11-97)

xb = x/b;
xr = x/r;
sqb = sqrt(xb.^2 +1);
sqr = sqrt(xr.^2 +1);

s = ( (xr.^2 + xb.^2 + log((abs(xb)+sqb)./(abs(xr)+sqr)) ).*sign(x) - ...
	xr.*sqr - xb.*( sqb + log((sqb-1)./(sqb+1)) ) )/(2*pi);