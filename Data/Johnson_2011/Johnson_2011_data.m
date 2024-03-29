function [ fig ] = Johnson_2011_data( plotflag )
%JOHNSON_2011_DATA 

if nargin < 1
    plotflag = false;
end

%% Figure 1A

% -------------------------------------------------------------------------
fig(1).A1 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'displacement', ...
    'yunit', 'um', ...
    'axis', {{}}, ...
    'description', '');

fig(1).A1.datasets(1) = struct(...
    'name', 'X', ...
    'data', csvread('Fig_1A_X.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

% -------------------------------------------------------------------------
fig(1).A2 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'Current', ...
    'yunit', 'nA', ...
    'axis', {{}}, ...
    'description', '');

fig(1).A2.datasets(1) = struct(...
    'name', 'I', ...
    'data', csvread('Fig_1A_I.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

%% Figure 1B

% -------------------------------------------------------------------------
fig(1).B1 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'displacement', ...
    'yunit', 'um', ...
    'axis', {{}}, ...
    'description', '');

fig(1).B1.datasets(1) = struct(...
    'name', 'X', ...
    'data', csvread('Fig_1B_X.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

% -------------------------------------------------------------------------
fig(1).B2 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'Current', ...
    'yunit', 'nA', ...
    'axis', {{}}, ...
    'description', '');

fig(1).B2.datasets(1) = struct(...
    'name', 'I', ...
    'data', csvread('Fig_1B_I.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

%% Figure 1C

% -------------------------------------------------------------------------
fig(1).C = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'voltage', ...
    'yunit', 'mV', ...
    'axis', {{}}, ...
    'description', '');

fig(1).C.datasets(1) = struct(...
    'name', 'X', ...
    'data', csvread('Fig_1C_V.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

%% Figure 1D

% -------------------------------------------------------------------------
fig(1).D = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'voltage', ...
    'yunit', 'mV', ...
    'axis', {{}}, ...
    'description', '');

fig(1).D.datasets(1) = struct(...
    'name', 'X', ...
    'data', csvread('Fig_1D_V.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

%% Figure 1E

% -------------------------------------------------------------------------
fig(1).E = struct( ...
    'xlabel', 'displacement', ...
    'xunit', 'um', ...
    'ylabel', 'current', ...
    'yunit', 'nA', ...
    'axis', {{}}, ...
    'description', '');

fig(1).E.datasets(1) = struct(...
    'name', 'I', ...
    'data', csvread('Fig_1E_I.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

%% Figure 1E

% -------------------------------------------------------------------------
fig(1).F = struct( ...
    'xlabel', 'displacement', ...
    'xunit', 'um', ...
    'ylabel', 'current', ...
    'yunit', 'nA', ...
    'axis', {{}}, ...
    'description', '');

fig(1).F.datasets(1) = struct(...
    'name', 'I', ...
    'data', csvread('Fig_1F_I.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});
    
%% Figure 2A

% -------------------------------------------------------------------------
fig(2).A1 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'displacement', ...
    'yunit', 'a.u.', ...
    'axis', {{}}, ...
    'description', '');

fig(2).A1.datasets(1) = struct(...
    'name', 'D', ...
    'data', csvread('Fig_2A_D.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

% -------------------------------------------------------------------------
fig(2).A2 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'Current', ...
    'yunit', 'nA', ...
    'axis', {{}}, ...
    'description', '');

fig(2).A2.datasets(1) = struct(...
    'name', 'I', ...
    'data', csvread('Fig_2A_I.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

fig(2).A2.datasets(2) = struct(...
    'name', 'DHS', ...
    'data', csvread('Fig_2A_I_DHS.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

% -------------------------------------------------------------------------
fig(2).A3 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'Voltage', ...
    'yunit', 'mV', ...
    'axis', {{}}, ...
    'description', '');

fig(2).A3.datasets(1) = struct(...
    'name', 'V', ...
    'data', csvread('Fig_2A_V.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

%% Figure 2B

fig(2).B1 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'displacement', ...
    'yunit', 'a.u.', ...
    'axis', {{}}, ...
    'description', '');

fig(2).B1.datasets(1) = struct(...
    'name', 'D', ...
    'data', csvread('Fig_2B_D.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

fig(2).B2 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'Current', ...
    'yunit', 'nA', ...
    'axis', {{}}, ...
    'description', '');

fig(2).B2.datasets(1)= struct(...
    'name', 'I', ...
    'data', csvread('Fig_2B_I.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

fig(2).B2.datasets(2) = struct(...
    'name', 'DHS', ...
    'data', csvread('Fig_2B_I_DHS.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

fig(2).B3 = struct( ...
    'xlabel', 'Time', ...
    'xunit', 'ms', ...
    'ylabel', 'Voltage', ...
    'yunit', 'mV', ...
    'axis', {{}}, ...
    'description', '');

fig(2).B3.datasets(1) = struct(...
    'name', 'V', ...
    'data', csvread('Fig_2B_V.csv'), ...
    'plotopt', ...
        {{'linestyle', '-', 'marker', 'none'}});

%% Figure 2C

fig(2).C = struct( ...
    'xlabel', 'Characteristic Frequency', ...
    'xunit', 'kHz', ...
    'ylabel', 'MT Current', ...
    'yunit', 'nA', ...
    'axis', {{'XScale', 'log'}}, ...
    'description', 'MT current (mean +- SEM) versus CF for three gerbil locations (filledsquares) and two rat locations (filled circles), number of measurements indicated by each point.');

% x, y, yerr
fig(2).C.datasets(1) = struct(...
    'name', 'gerbil', ...
    'data', [ ...
        0.35, 1.18681, 0.11209
        0.9, 1.49670, 0.16814
        2.5, 2.18901, 0.08242 ], ...
    'plotopt', ...
       {{'linestyle', 'none', 'marker', 's'}});

% x, y, yerr
fig(2).C.datasets(2) = struct(...
    'name', 'rat', ...
    'data', [ ...
        4.0, 2.17582, 0.13847
        10.0, 2.91099, 0.13846 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 'o'}});

% x, y
fig(2).C.datasets(3) = struct(...
    'name', 'fit', ...
    'data', [ ...
        0.24909, 0.97305
       12.29745, 2.95666 ], ...
    'plotopt', ...
        {{'linestyle', '--', 'marker', 'none'}});


%% Figure 2D

fig(2).D = struct( ...
    'xlabel', 'Characteristic Frequency', ...
    'xunit', 'kHz', ...
    'ylabel', 'Resting Popen', ...
    'yunit', '', ...
    'axis', {{'XScale', 'log'}}, ...
    'description', 'Resting open probability (P open; mean +- SEM) of the MT channels for gerbil and rat locations (P7–P10 animals). Currents recorded at 84 mV holding potential, T = 22-24 deg C.');

% x, y, yerr
fig(2).D.datasets(1) = struct(...
    'name', 'gerbil', ...
    'data', [ ...
        0.35, 0.42825, 0.03992
        0.9, 0.45093, 0.04622
        2.5, 0.58913, 0.07774 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 's'}});


% x, y, yerr
fig(2).D.datasets(2) = struct(...
    'name', 'rat', ...
    'data', [ ...
        4.0, 0.42924, 0.04201
       10.0, 0.38680, 0.05462 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 'o'}});


% x, y
fig(2).D.datasets(3) = struct(...
    'name', 'fit', ...
    'data', [ ...
        0.20703, 0.45790
       17.26792, 0.45588 ], ...
    'plotopt', ...
        {{'linestyle', '--', 'marker', 'none'}});

    
    
%% Figure 6C

fig(6).C = struct( ...
    'xlabel', 'Characteristic Frequency', ...
    'xunit', 'kHz', ...
    'ylabel', 'Resting potential', ...
    'yunit', 'mV', ...
    'axis', {{'XScale', 'log'}}, ...
    'description', '');

% x, y, yerr
fig(6).C.datasets(1) = struct(...
    'name', 'gerbil', ...
    'data', [ ...
         350, -26.980, 0.8580
         900, -38.050, 2.9080
        2500, -43.019, 2.9560 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 's'}});


% x, y, yerr
fig(6).C.datasets(2) = struct(...
    'name', 'rat', ...
    'data', [ ...
         4000, -42.071, 2.8610
        10000, -50.376, 1.8120 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 'o'}});


% x, y
fig(6).C.datasets(3) = struct(...
    'name', 'fit', ...
    'data', [ ...
        207.62, -40.038
        17149, -40.083 ], ...
    'plotopt', ...
        {{'linestyle', '--', 'marker', 'none'}});

    
%% Figure 6D

fig(6).D = struct( ...
    'xlabel', 'Characteristic Frequency', ...
    'xunit', 'kHz', ...
    'ylabel', 'Conductance', ...
    'yunit', 'nS', ...
    'axis', {{'XScale', 'log'}}, ...
    'description', '');

% x, y, yerr
fig(6).D.datasets(1) = struct(...
    'name', 'gerbil', ...
    'data', [ ...
         350, 26.919, 3.2430
         900, 51.081, 3.8920
        2500, 87.730, 9.7290 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 's'}});


% x, y, yerr
fig(6).D.datasets(2) = struct(...
    'name', 'rat', ...
    'data', [ ...
         4000, 79.946, 8.9190
        10000, 149.84, 12.160 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 'o'}});


% x, y
fig(6).D.datasets(3) = struct(...
    'name', 'fit', ...
    'data', [ ...
        322.10, 22.703
        14455, 150.65 ], ...
    'plotopt', ...
        {{'linestyle', '--', 'marker', 'none'}});
    

%% Figure 7A

fig(7).A = struct( ...
    'xlabel', 'Characteristic Frequency', ...
    'xunit', 'Hz', ...
    'ylabel', 'Capacitance', ...
    'yunit', 'pF', ...
    'axis', {{'XScale', 'log'}}, ...
    'description', '');

% x, y, yerr
fig(7).A.datasets(1) = struct(...
    'name', 'gerbil', ...
    'data', [ ...
         350, 18.34936, 1.5109
         900, 17.44937, 1.8847
        2500, 13.85557, 0.5452 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 's'}});


% x, y, yerr
fig(7).A.datasets(2) = struct(...
    'name', 'rat', ...
    'data', [ ...
         4000, 11.00007,  1.3084
        10000, 5.30360, 1.0125 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 'o'}});

%% Suplementary Figure 1A

fig(9).A = struct( ...
    'xlabel', 'Characteristic Frequency', ...
    'xunit', 'Hz', ...
    'ylabel', 'Conductance', ...
    'yunit', 'nS', ...
    'axis', {{'XScale', 'log'}}, ...
    'description', '');

% x, y
fig(9).A.datasets(1) = struct(...
    'name', 'gerbil', ...
    'data', [ ...
          250, 28.99313
          900, 56.76732
         2500, 90.31039
        16000, 255.97451 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 's'}});

% x, y
fig(9).A.datasets(2) = struct(...
    'name', 'rat', ...
    'data', [ ...
         4000, 81.55966
        10000, 240.97235 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 'o'}});

%% Suplementary Figure 1B

fig(9).B = struct( ...
    'xlabel', 'Characteristic Frequency', ...
    'xunit', 'Hz', ...
    'ylabel', 'Conductance', ...
    'yunit', 'nS', ...
    'axis', {{'XScale', 'log'}}, ...
    'description', '');

% x, y, yerrmin, yerrmax
fig(9).B.datasets(1) = struct(...
    'name', 'gerbil', ...
    'data', [ ...
          350, 22.63023, 2.2943, 2.2854
          900, 28.50805, 3.1626, 3.2163
         2500, 41.62104, 2.4184, 2.4733
        16000, 88.52689, 0, 0 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 's'}});

% x, y
fig(9).B.datasets(2) = struct(...
    'name', 'rat', ...
    'data', [ ...
          4000, 41.62104, 3.0008, 2.4733
         10000, 55.31005, 2.9906, 2.9121
         20000, 65.90851, 3.4303, 3.3219 ], ...
    'plotopt', ...
        {{'linestyle', 'none', 'marker', 'o'}});

    
%%

if plotflag
    
    for i = 1:numel(fig)

        s = fields(fig(i));

        for j = 1:numel(s)

            if ~isempty(fig(i).(s{j}))

                figure
                hold on
    
                Johnson_2011_plot( fig(i), s{j} )
                
                title(sprintf('Johnson 2011, Fig. %d, %s', i, s{j}))
                
            end
        end
    end
    
end

end

