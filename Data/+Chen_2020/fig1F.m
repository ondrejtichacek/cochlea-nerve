function [x, y, eplus, eminus] = fig1F()
arguments
end

folder = fileparts(mfilename('fullpath'));

f = fullfile(folder, 'fig1F_noise.csv');
tmp = readmatrix(f);

x = tmp(:,1);
y = tmp(:,2);

f = fullfile(folder, 'fig1F_noise_err.csv');
tmp = readmatrix(f);
eplus = -y + tmp(1:2:end,2);
eminus = y - tmp(2:2:end,2);


end

