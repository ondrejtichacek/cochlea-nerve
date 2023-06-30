function [x, y] = fig2b()
%FIG2 

folder = fileparts(mfilename('fullpath'));

f2b = fullfile(folder, 'f2b.csv');

tmp = readmatrix(f2b);
x = tmp(:,1);
y = tmp(:,2);

end

