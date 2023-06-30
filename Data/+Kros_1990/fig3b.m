function [x, y, xl, yl] = fig3b()
arguments
end

folder = fileparts(mfilename('fullpath'));

f_data = fullfile(folder, 'f3b1.csv');
f_fit = fullfile(folder, 'f3b2.csv');

tmp = readmatrix(f_data);
x = tmp(:,1);
y = tmp(:,2);


tmp = readmatrix(f_fit);
xl = tmp(:,1);
yl = tmp(:,2);


end

