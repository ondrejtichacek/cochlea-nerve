function [x, y] = fig2c()
%FIG2 

folder = fileparts(mfilename('fullpath'));

f2c = fullfile(folder, 'f2c.csv');

tmp = readmatrix(f2c);
x = tmp(:,1);
y = tmp(:,2);

end

