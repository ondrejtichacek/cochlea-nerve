% Startup script setting up the path for CochleNerve model

% Set UTF-8 character encoding
addpath(fullfile('ext', 'IsOctave'));

if IsOctave
else
    feature('DefaultCharacterSet', 'UTF8');
end

% Restore default path
% restoredefaultpath;

% Add the genpath_exclude function to path
addpath(fullfile('externals', 'genpath_exclude'));

try
    % Generate and add path excluding the directories '.git' and '#.*'
    % and directories containing a file named #UNUSED

    addpath(genpath_exclude(pwd, {'\.git','#.*'}, {'#UNUSED'}));
    
catch ME
    disp(ME)
    warning('genpath_exclude error - switching to original genpath!')
    addpath(genpath(pwd));
end


if usejava('jvm') && ~feature('ShowFigureWindows')
    % text mode
else
    % GUI mode
    set(0, 'DefaultFigureWindowStyle', 'docked')
    set(groot, 'defaultLineLineWidth', 1)
end