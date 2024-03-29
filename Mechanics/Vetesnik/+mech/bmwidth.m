function [ b, bmwdata, hfig ] = bmwidth( x, plotflag )
%BMWIDTH
% x ... normalized coordinate
% b ... BM width
arguments
    x
    plotflag (1,1) logical = false
end

% variable explanation:
%   bmwdata(:,1) = abscissae,
%   bmwdata(:,2) = ordinates

bmwdata = [ ...
   0.0000000e+00   1.7991004e-03
   9.4696969e-03   2.4756852e-03
   1.5151515e-02   2.8139776e-03
   2.8409091e-02   3.3214162e-03
   4.1666666e-02   3.5582210e-03
   5.8712120e-02   3.8288549e-03
   8.3333334e-02   4.1107652e-03
   2.2916667e-01   5.4921256e-03
   3.2954546e-01   6.3860477e-03
   3.8446969e-01   6.9580594e-03
   4.3560606e-01   7.4936891e-03
   5.7954546e-01   9.0723868e-03
   6.8560606e-01   1.0425557e-02
   7.3863637e-01   1.1271288e-02
   7.9545454e-01   1.2398929e-02
   8.3049241e-01   1.3301042e-02
   8.6553031e-01   1.4372301e-02
   8.8068183e-01   1.4710593e-02
   8.9583334e-01   1.4879740e-02
   9.0909091e-01   1.4950004e-02
   9.2424243e-01   1.4936122e-02
   9.4318183e-01   1.4654211e-02
   9.5833334e-01   1.4275919e-02
   9.7916666e-01   1.3413806e-02
   9.8484849e-01   1.3160087e-02
   9.9431817e-01   1.2511693e-02
   1.0000000e+00   1.2060636e-02
];

% b = interpol(bmwdata(:,1),bmwdata(:,2),x);
F = griddedInterpolant(bmwdata(:,1),bmwdata(:,2));
b = F(x);

if plotflag
    
    BETA = 33.5; % [mm]
    
    hfig = figure();
    
	plot(bmwdata(:,1), BETA*bmwdata(:,2), bmwdata(:,1), BETA*bmwdata(:,2), ':r'),
    
    %axis([0 1 0 0.02*BETA]);
    
    title('Basilar membrane width');
	xlabel('Fractional distance from stapes');
	ylabel('mm');
    
end

end