function [val] = canplot()
%CANPLOT 

val = ~(usejava('jvm') && ~feature('ShowFigureWindows'));

end

