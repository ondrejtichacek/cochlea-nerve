function handles = parfeval_ui_monitor(pf, identifiers, args)
arguments
    pf
    identifiers
    args.handles = []
end
    
% f = uifigure();
% tabgp = uitabgroup(f);
% t = uitab(tabgp);

if isempty(args.handles)

    f = figure();

    drawnow()

    uit_1 = uitable(f);
    uit_1.Position = [20, 40, f.Position(3) - 40, f.Position(4)/2 - 40];

    uit_2 = uitable(f);
    uit_2.Position = [20, 20 + f.Position(4)/2, f.Position(3) - 40, f.Position(4)/2 - 40];
    
    hb = uicontrol( f, ...
        'Style', 'PushButton', ...
        'String', 'Break batchGathering:fetchNext', ...
        'Callback', 'delete(gcbf)');
    
    handles = struct( ...
        'figure', f, ...
        'stop_button', hb, ...
        'uitable_1', uit_1, ...
        'uitable_2', uit_2 ...
        );
    
else
    handles = args.handles;
end

%%

vars = {};
vars = [vars, {pf.ID}'];
vars = [vars, {pf.State}'];
vars = [vars, cellfun(@(x) formatError(x), {pf.Error}', 'UniformOutput', false)];

var_names = {'ID', 'State', 'Error'}; %, 'CreateDateTime', 'StartDateTime', 'FinishDateTime'};

id_names = fieldnames(identifiers);

id_names = id_names(1:min([2, numel(id_names)]));

for i = 1:numel(id_names)
    vars = [vars, {identifiers.(id_names{i})}'];
end

var_names = [var_names, id_names'];

vars = [vars, cellfun(@(x) formatDiary(x), {pf.Diary}', 'UniformOutput', false)];

var_names = [var_names, {'Diary'}];

running = strcmp(vars(:,2),'running');

%%

data = vars;
data = data(running,:);

for i = 1:numel(data)
    if isa(data{i}, 'cell')
        data{i} = 'can''t display cell';
    end
        
end

handles.uitable_1.Data = data;

handles.uitable_1.ColumnName = var_names;
handles.uitable_1.ColumnFormat = ([{[] [] []}, repmat({[]}, 1, numel(id_names)), {'char'}]);
% handles.uitable_1.ColumnWidth = 'fit';
handles.uitable_1.ColumnWidth = ([{40 80 80}, repmat({60}, 1, numel(id_names)), {720}]);


%%

data = vars;    
data = data(~running,:);

for i = 1:numel(data)
    if isa(data{i}, 'cell')
        data{i} = 'can''t display cell';
    end
        
end

handles.uitable_2.Data = data;

handles.uitable_2.ColumnName = var_names;
handles.uitable_2.ColumnFormat = ([{[] [] []}, repmat({[]}, 1, numel(id_names)), {'char'}]);
% handles.uitable_2.ColumnWidth = 'fit';
handles.uitable_2.ColumnWidth = ([{40 80 80}, repmat({60}, 1, numel(id_names)), {720}]);

%%

% s = uistyle('BackgroundColor',[1 0.6 0.6]);
% addStyle(uit,s,'cell',[row, col]);
    
end

function [ line ] = formatError(E)
if isempty(E)
    line = '';
else
    if isa(E, 'MException')
        line = sprintf('%s: %s', E.identifier, E.message);
    else
        line = sprintf('WARN: class(E) = %s', class(E));
    end
end
end

function [ line ] = formatDiary(str)
s = regexp(str,'\n','split');
if numel(s) >= 2
    line = s{end-1};
    line = regexprep(line,'\s',' ');
    line = regexprep(line,'→',' ');
else
    line = '';
end
end