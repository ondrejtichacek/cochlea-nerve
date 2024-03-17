function [ xp, yp, xm, ym, ap, am ] = cochlearduct()
%COCHLEARDUCT

% disp('	Reading human cochlea data...');

load +mech/cddata.mat xp yp xm ym ap am

% File contains cochlear duct profile data
%   All lengths are scaled to 1 unit == 33.5 mm
%   Data from the upper side:  xp, yp
%   Data from the lower side:  xm, ym
%   ap, am = cross--sectional areas
%   "p" is for upper profile, "m" is for lower.

end