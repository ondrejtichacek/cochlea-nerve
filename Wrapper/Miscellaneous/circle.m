function [xp, yp] = circle(r)
%CIRCLE 
ang = linspace(0, 2*pi); 
xp = r*cos(ang);
yp = r*sin(ang);
end