function err=intgausserr(X0,x,y)

thrsh	 = X0(1);
noisestd = X0(2);

f = 1/2*(erf((x-thrsh)/sqrt(2)./noisestd)+1);

err = norm(f-y);