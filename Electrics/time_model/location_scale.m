function [loc_scale] = location_scale(F, A, gain, xgrid, width, args)
%LOCATION_SCALE 
arguments
    F
    A
    gain
    xgrid
    width
    
    args.plotflag (1,1) logical = false
end

[cx, cf] = characteristic_frequency_model('Li');
x2 = cx(F);
[cx, cf] = characteristic_frequency_model('Greenwood');
x1 = cx(F);

% Equivalent amplitude takes into account the gain, which sharpens the
% response. Add up to 80 dB if the amplifier is off

Aeq = A + 80 * (1 - gain);

A1 = 30;
A2 = 80;

if Aeq < A1
    xc = x2;
    
    w = [width, width];

elseif Aeq > A2
    xc = x1;
    w = [width + 0.1, width];

else
    xc = x2 - (x2 - x1) * (Aeq - A1) / (A2 - A1);
    w = [width + 0.1 * (Aeq - A1) / (A2 - A1), width];
end

a = 0.9;

loc_scale = 1 - a  ...
    + gaussian(xgrid, a, xc, w(1)) .* (xgrid <= xc) ...
    + gaussian(xgrid, a, xc, w(2)) .* (xgrid > xc);

loc_scale(loc_scale > 1) = 1;

if args.plotflag
    figure
    plot(xgrid, loc_scale)
end

end

