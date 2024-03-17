function [m, hfig] = mass(x, plotflag)
%MASS
%  Computes the effective organ of Corti mass per unit BM length as a function of
%  fractional distance from stapes X. This function loads and interpolates BM-width
%  data from BMWDATA.MAT. The unit dimension of M is Kg/BETA, where
%  BETA = 33.5 mm is the BM length. If PLOTFLAG = 1 the mass density profile
%  is plotted [def. PLOTFLAG = 0]
%

if nargin < 2
    plotflag = 0;
end

if nargin < 1
    x = linspace(0,1,200);
    plotflag = 1;
end

x = x(:);

%% UNIT CONVERSION

BETA = 33.5; % [mm]
METER = 1e3 / BETA; % 1 m in BETA
micron = METER*1e-6; % 1 micron in BETA

% rho = 0.0376; % = 1000/(29.85)^3, water density in Kg/BETA^3, as 1m = 29.85 BETA
rho = 1000/METER^3;

%%

b = mech.bmwidth(x);

% On account of radial BM curvature, the effective BM width is assumed geometric-width/2
b = b/2;

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ho = 0.02/BETA;    %  Mean organ of Corti  thickness in BETA units (it is assumed constant)
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

h0 = 20*micron; % 20 micron
% putative height of the organ of Corti at base
% (effective for shear viscosity)

h = h0*exp(log(4)*x);
% effective organ of Corti height
% (increases by 4 times from base to apex)

w0 = 50*micron; % 50 micron,
% putative width of the organ of Corti at base
% (effective for shear viscosity)

w = w0*exp(log(4)*x);
% effective organ of Corti width
% (increases by 4 times from base to apex)

w = sqrt(w.*b);    % geometrical mean of effective widths

m = rho*(h.*b);    % vector  m  is a row
m = (m + [m(1); m(1:end-1)])/2;  % Local averaging


if plotflag == 1
    hfig = figure();
    plot(x,m);
    title('Organ of Corti mass per unit BM length');
    xlabel('Fractional distance from stapes');
    ylabel('Kg/BETA (BETA = 33.5 mm).');
end

end