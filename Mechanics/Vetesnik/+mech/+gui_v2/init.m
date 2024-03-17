function [ mechguiopt ] = init( x, t, Sig, mnaopt, args )
%INIT
arguments
    x
    t
    Sig
    mnaopt (1,1) mnaOpt
    args.frequency = []
    args.init_state (1,:) char = 'run'
    args.t0 = []
    args.tf = []
end

if isempty(args.t0)
    args.t0 = t(1);
end
if isempty(args.tf)
    args.tf = t(end);
end


draw_patches = false;
draw_patches = true;

N = numel(x);

NumberOfSubfigures = 5;
NumberOfSubfigures = 7;
NumberOfSubfigures = 8;
NumberOfSubfigures = 10;
NumberOfSubfigures = 12;

time_unit = 'ms';

if ~draw_patches
    p = [];
end

%% figure
hfig = figure( ...
    'Name','TW model', ...
    'DoubleBuffer','on', ...
    'NumberTitle','off'); % figure handle

set(hfig, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.2, 0.6, 0.6]);

sp_cnt = 0;

%% signal
sp_size = 1;
ax_sig = subplot(NumberOfSubfigures,1,sp_cnt + (1:sp_size));
sp_cnt = sp_cnt + sp_size;

pos = ax_sig.Position;
pos(3) = 0.75;

set(ax_sig, 'Xlim', [args.t0.(time_unit), args.tf.(time_unit)], ...            %'Ylim', [-Ma1, Ma1], ...
            'Box', 'on', ...
            'Position', pos);

xlabel(sprintf('Time (%s)', time_unit));
ylabel('Amplitude');

c = [0.7 0.7 0.7
     0   0   1];

for ii = 1:size(Sig,1)
    line_sig(ii) = line( ... % signal plot
        'XData', t.(time_unit), ...
        'YData', Sig(ii,:), ...
        'Color', c(ii,:), ...
        'Parent', ax_sig);
end

% Ma1 = -AM*om2;
line_0 = line( ... % vertical line -- current position
    'Xdata', [args.t0.(time_unit), args.t0.(time_unit)], ...
    'Ydata', ax_sig.YLim, ... %     'Ydata', [-Ma1, Ma1], ...
    'Color', 'r', ...
    'Parent',ax_sig);

%% spectrogram
sp_size = 1;
ax_spectrogram = subplot(NumberOfSubfigures,1,sp_cnt + (1:sp_size));
sp_cnt = sp_cnt + sp_size;

pos = ax_spectrogram.Position;
pos(3) = 0.75;

xlabel(sprintf('Time [%s]', time_unit));
ylabel('Amplitude');

fs = 1/mean(diff(t));

signal = Sig(1,:);

% window = min(t(end)/2, Time(20, 'ms')) * fs;
window = round(numel(signal)/4.5/10);

if numel(signal) > 2
    spectrogram(signal,window,[],[],fs.Hz,'yaxis');
    ylim([0,20])
    colorbar off
end

set(ax_spectrogram, 'Xlim', [args.t0.(time_unit), args.tf.(time_unit)], ...            %'Ylim', [-Ma1, Ma1], ...
            'Box', 'on', ...
            'Position', pos);

% c = [0.7 0.7 0.7
%      0   0   1];
% 
% for ii = 1:size(Sig,1)
%     line_sig(ii) = line( ... % signal plot
%         'XData', t.(time_unit), ...
%         'YData', Sig(ii,:), ...
%         'Color', c(ii,:), ...
%         'Parent', ax_sig);
% end
% 
% % Ma1 = -AM*om2;
% line_0 = line( ... % vertical line -- current position
%     'Xdata', [t(1).(time_unit), t(1).(time_unit)], ...
%     'Ydata', ax_sig.YLim, ... %     'Ydata', [-Ma1, Ma1], ...
%     'Color', 'r', ...
%     'Parent',ax_sig);

%% Heatmap 1

sp_size = 3;

for kk = 3

    ax(kk) = subplot(NumberOfSubfigures,1,sp_cnt + (1:sp_size));
    sp_cnt = sp_cnt + sp_size;

    set(ax(kk), ...
        'Color', [0.9,0.9,0.9], ...
        'YDir', 'reverse', ...
        'Box', 'on', ...
        'Ylim', [0, 1], ...
        'Xlim', [args.t0.(time_unit), args.tf.(time_unit)]);
    yticks([])
    ylabel('rel. pos.')
    xlabel('Time (ms)')

    yy = mnaopt.xgrid;
    xx = t.(time_unit);

    V = zeros(numel(xx), numel(yy));

    hold on

    imagesc('XData', xx, 'YData', yy, 'CData', V)
    % [XX, YY] = meshgrid(xx, yy);
    % pcolor(XX', YY', V)
    
    colormap(ax(kk), flip(gray))
    % colormap(ax(kk), seismic)

    pos = ax(kk).Position;
    pos(3) = 0.75;
    ax(kk).Position = pos;

    plot(-1,-1, 'r.');

    plot(-1,-1, 'k.');

    text(0.7,0.9, 'IHC receptor pot.', Units='normalized', FontSize=14, Color='k')

    % legend('Interpreter', 'none') % legend is super slow    

%     pos = ax(kk).Position;
%     pos(3) = 0.75;
%     ax(kk).Position = pos;
end

%%  Heatmap 2

sp_size = 3;

for kk = 4

    ax(kk) = subplot(NumberOfSubfigures,1,sp_cnt + (1:sp_size));
    sp_cnt = sp_cnt + sp_size;

    set(ax(kk), ...
        ...'Color', [1,1,1], ...
        'Color', [0.9,0.9,0.9], ...
        'YDir', 'reverse', ...
        'Box', 'on', ...
        'Ylim', [0, 1], ...
        'Xlim', [args.t0.(time_unit), args.tf.(time_unit)]);
    yticks([])
    ylabel('rel. pos.')
    xlabel('Time (ms)')

    yy = mnaopt.xgrid;
    xx = t.(time_unit);

    V = zeros(numel(xx), numel(yy));

    hold on

    imagesc('XData', xx, 'YData', yy, 'CData', V)
    % [XX, YY] = meshgrid(xx, yy);
    % pcolor(XX', YY', V)
    
    cmap = hot; cmap = cmap(1:end-20,:);
    cmap = [cmap; [1,1,1]];
    colormap(ax(kk), flip(cmap))
    % colormap(ax(kk), seismic)

    pos = ax(kk).Position;
    pos(3) = 0.75;
    ax(kk).Position = pos;

    plot(-1,-1, 'k.');

    text(0.7,0.9, 'Neural activity', Units='normalized', FontSize=14, Color='k')
end

%% displacement

[cx_1, ~] = characteristic_frequency_model('Greenwood');
[cx_2, ~] = characteristic_frequency_model('Li');

sp_size = 2;

for kk = 1:2
    ax(kk) = subplot(NumberOfSubfigures,1,sp_cnt + (1:sp_size));
    % legend('Interpreter', 'none') % legend is super slow
    sp_cnt = sp_cnt + sp_size;

    pos = ax(kk).Position;
    pos(3) = 0.75;
    ax(kk).Position = pos;
    
    
    if ~isempty(args.frequency)
        cx = cx_1(args.frequency.Hz);
        xline(cx)
        cx = cx_2(args.frequency.Hz);
        xline(cx)
    end

    NumLines = 2;
    % line_colors = lines(NumLines);
    line_colors_r = [1,0,0; 1,0,0.5];
    line_colors_b = [0,0,1; 0.5,0,1];

    for i = 1:NumLines

        yyaxis left

        l(kk,i,1) = line( ...
            'XData', x, ...
            'YData', zeros(N,1), ...
            'Color', line_colors_b(i,:), ...
            'Parent', ax(kk));
        
        text_label(kk,i,1) = text( ...
            0.1, 0.9-0.8*(i-1), '', ...
            Units='normalized', ...
            FontSize=14, ...
            Color=line_colors_b(i,:));

        if draw_patches
            lE = NaN*-eps*ones(1,N);
            uE = NaN*eps*ones(1,N);
        
            yP = [lE, fliplr(uE)];
            xP = [x, fliplr(x)];
        
            p(kk,i,1) = patch(xP,yP,1);
        
            set(p(kk,i,1), ...
                'facecolor',line_colors_b(i,:), ...
                'edgecolor','none', ...
                'facealpha',0.1)
        end

        set(ax(kk), ...
            'XLim',[x(1),x(N)], ...
            'YLim', [-1e3*eps, 1e3*eps], ...
            'YLimMode', 'manual', ...
            'Box', 'off');

        % set(ax(kk), 'yscale', 'log')
        
        yyaxis right

        l(kk,i,2) = line( ...
            'XData', x, ...
            'YData', zeros(N,1), ...
            'Color', line_colors_r(i,:), ...
            'Parent', ax(kk));
        
        text_label(kk,i,2) = text( ...
            0.7,0.9-0.8*(i-1), '', ...
            Units='normalized', ...
            FontSize=14, ...
            Color=line_colors_r(i,:));

        if draw_patches
            lE = NaN*-eps*ones(1,N);
            uE = NaN*eps*ones(1,N);
            
            yP=[lE, fliplr(uE)];
            xP=[x, fliplr(x)];
            
            p(kk,i,2) = patch(xP,yP,1);
            
            set(p(kk,i,2), ...
                'facecolor',line_colors_r(i,:), ...
                'edgecolor','none', ...
                'facealpha',0.1)
        end

        set(ax(kk), ...
            'XLim',[x(1),x(N)], ...
            'YLim', [-1e3*eps, 1e3*eps], ...
            'YLimMode', 'manual', ...
            'Box', 'off');
        
        xlabel('Normalized distance from stapes')

        %set(ax(kk), 'yscale', 'log')
        
    end

    ax(kk).YAxis(1).Color = [0 0 1];
    ax(kk).YAxis(2).Color = [1 0 0];

end

text_label(1,1,1)% = text(0.7,0.9, 'BM displacement', Units='normalized', FontSize=14, Color='b');
text_label(2,1,2)
%text_label = text(0.7,0.9, 'IHC receptor potential', Units='normalized', FontSize=14, Color='r')

%%

axes(ax(3)) % for some reason, it would not update imagesc otherwise

%% twa
%Ma = 1.1*max(TWa);
% Ma = eps;
% TWa = zeros(N,1);
% ax_twa = axes('Position', ax(1).Position, ...
%             'XLim',  [x(1),x(N)], ...
%             'YLim',[-Ma,Ma], ...
%             'Box', 'off', ...
%             'XAxisLocation','top',...
%             'YAxisLocation','right', ...
%             'XTick', [], ...
%             'Color','none');
%         
% line_twa = line('XData', x, ...
%             'YData', TWa, ...
%             'Color', 0.5*[1,1,1], ...
%             'Parent', ax_twa);

%%

mechguiopt = mechGuiOpt( ...
    'Numstacks', N, ...
    'last_t', 0, ...
    'H_FIG', hfig, ...
    'ax_sig', ax_sig, ...
    'ax', ax, ...
    ...'ax_twa', ax_twa, ...
    'line_sig', line_sig, ....
    'line_0', line_0, ...
    'lines', l, ...
    'text_label', text_label, ...
    'draw_patches', draw_patches, ...
    'patches', p ...
    );
    ...'line_twa', line_twa);

%% GUI interaction

channels = {mnaopt.channels.name};

signals = { ...
    'select', ...
    'F_OHC', ...
    'mech_BMx', 'mech_TMx', ...
    'mech_BMv', 'mech_TMv', ...
    'stapes_mech_BMx', 'stapes_mech_TMx', 'stapes_mech_BMv', 'stapes_mech_TMv', ...
    'amplifier_mech_BMv', 'amplifier_mech_TMv', ...
    'amplification_BM', 'amplification_TM', ...
    'curr_IHC', 'curr_OHC', ...
    'volt_IHC', 'volt_OHC', ...
    'volt_IHC_ic', 'volt_IHC_ec', 'volt_OHC_ic', 'volt_OHC_ec', ......
    'tol_mech_BMx', 'tol_mech_TMx', 'tol_mech_BMv', 'tol_mech_TMv', ...
    'tol_curr_IHC', 'tol_curr_OHC', ...
    'tol_volt_IHC_ic', 'tol_volt_IHC_ec', 'tol_volt_OHC_ic', 'tol_volt_OHC_ec', ......
    };

for i = 1:numel(channels)
    signals{end+1} = ['p_channel_', channels{i}];
    signals{end+1} = ['dp_channel_', channels{i}];
end

resistors = mnaopt.circuit.resistors;
for ii = 1:numel(resistors)
    if resistors(ii).special == "radial"
        signals{end+1} = ['R_', resistors(ii).name];
    end
end

def_signals = {
    'mech_BMx', 'select', ...'amplification_BM', ...
    'select', 'volt_IHC', ...
    ...'mech_BMx', 'mech_TMx', ...
    ...'mech_BMv', 'volt_OHC', ...
    ...'amplifier_mech_BMv', 'amplification_BM', ...
    };

def_modifiers = {
    'Absolute', 'Absolute', ...
    'Absolute', 'Relative', ...
    };

fig = mechguiopt.H_FIG;

heat_signals = { ...
    'select', ...
    'V', ...
    'synchronization', ...
    };

def_heat_signals = {
    'V', ...
    'synchronization', ...
    };

heat_over_signals = { ...
    'select', ...
    'release_event', ...
    'action_potential', ...
    };

def_heat_over_signals = { ...
    'select', ...'release_event', ...
    'action_potential', ...
    };

% -------------------------------------------------------------------------

jj = 1;

for kk = 1:2
    for ll = 1:1 % left/right
        for ii = 1:1 % num signals
        
        pos = [0.005 + 0.90*(ll-1) 0.7 - 0.10*(ii-1) - 0.10*(kk-1) 0.09 0.05];

        heat_drop(kk,ii,ll) = uicontrol(fig, ...
            'Style', 'popupmenu', ...
            'Units', 'normalized', ...
            'Position', pos, ...
            'String', heat_signals);

        pos = pos - [0 0.05 0 0];
        
        heat_over_drop(kk,ii,ll) = uicontrol(fig, ...
            'Style', 'popupmenu', ...
            'Units', 'normalized', ...
            'Position', pos, ...
            'String', heat_over_signals);

        end
        
        ind = find(strcmp(def_heat_signals{jj}, heat_signals), 1);
        heat_drop(kk,1,ll).Value = ind;

        ind = find(strcmp(def_heat_over_signals{jj}, heat_over_signals), 1);
        heat_over_drop(kk,1,ll).Value = ind;

        jj = jj + 1;

    end
end

% -------------------------------------------------------------------------

jj = 1;

for kk = 1:2
    for ll = 1:2 % left/right
        for ii = 1:2
        
        pos = [0.005 + 0.90*(ll-1) 0.45 - 0.10*(ii-1) - 0.25*(kk-1) 0.09 0.05];

        sig_drop(kk,ii,ll) = uicontrol(fig, ...
            'Style', 'popupmenu', ...
            'Units', 'normalized', ...
            'Position', pos, ...
            'String', signals);

        pos = pos - [0 0.05 0 0];
        
        rel_drop(kk,ii,ll) = uicontrol(fig, ...
            'Style', 'popupmenu', ...
            'Units', 'normalized', ...
            'Position', pos, ...
            'String', ["Absolute", "Relative"]);

            if ii == 2
    
                pos = pos - [0 0.05 0 0];
        
                log_scale(kk,ii,ll) = uicontrol(fig, ...
                    'Style', 'checkbox', ...
                    'Units', 'normalized', ...
                    'Callback', {@log_scale_callback, ax(kk), ll}, ...
                    'Position', pos, ...
                    'String', {'Log'});
            end
        end
        
        ind = find(strcmp(def_signals{jj}, signals), 1);
        % ind = jj + 1;
        sig_drop(kk,1,ll).Value = ind;

        ind = find(strcmp(def_modifiers{jj}, ["Absolute", "Relative"]), 1);
        % ind = jj + 1;
        rel_drop(kk,1,ll).Value = ind;
        
        jj = jj + 1;

    end
end


% -------------------------------------------------------------------------

switch args.init_state
    case 'run'
        state_str = 'Pause';
    case 'paused'
        state_str = 'Resume';
    otherwise
        error('Unknown option %s', args.init_state)
end

set(ax_sig, 'ButtonDownFcn', {@mouse_click_signal_callback, mechguiopt})

btn(1) = uicontrol(fig, ...
    'Style','pushbutton', ...
    'Units', 'normalized', ...
    'Callback', {@pause_resume_callback, mechguiopt}, ...
    'Position', [0.005 0.92 0.08 0.05], ...
    'String', state_str);

btn(2) = uicontrol(fig, ...
    'Style','pushbutton', ...
    'Units', 'normalized', ...
    'Callback', {@reset_scale_callback, mechguiopt}, ...
    'Position', [0.005 0.85 0.08 0.05], ...
    'String', 'Rescale');

btn(3) = uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Callback', {@abort_callback, mechguiopt}, ...
    'Position', [0.92-0.005 0 0.08 0.05], ...
    'String', {'Abort'});

c = uicontrol(fig, ...
    'Style', 'popupmenu', ...
    'Units', 'normalized', ...
    'Callback', {@update_frequency_callback, mechguiopt}, ...
    'Position', [0.92-0.005 0.92 0.08 0.05], ...
    'String', {'100 kHz', '50 kHz', '20 kHz', '10 kHz', '1 kHz', 'inf',});

c.Value = 2; % 50 kHz
update_frequency_callback(c, [], mechguiopt); % force initial value

% set(mechguiopt.H_FIG, 'WindowKeyPressFcn', @KeyPress);
mechguiopt.H_FIG.UserData = mechguiopt;

mechguiopt.uicontrol = struct( ...
    'heat_drop', heat_drop, ...
    'heat_over_drop', heat_over_drop, ...
    'signal', sig_drop, ...
    'rel_switch', rel_drop, ...
    'btn_start', btn(1), ...
    'btn_reset', btn(2), ...
    'btn_abort', btn(3), ...
    'update_freq', c ...
    );
end


function mouse_click_signal_callback(src, eventdata, mechguiopt)
    coordinates = src.CurrentPoint(1,1:2);

    mechguiopt.loader_t_click = coordinates(1) / 1e3; % in s
end

function log_scale_callback(hObject, eventdata, ax, ll)

    if ll == 1
        yyaxis(ax, 'left')
    elseif ll == 2
        yyaxis(ax, 'right')
    end
    if hObject.Value == 1
        set(ax, 'yscale', 'log')
    else
        set(ax, 'yscale', 'linear')
    end

end

function pause_resume_callback(hObject, eventdata, handles)

    button_state = hObject.String;

    switch button_state
        case 'Pause'
            hObject.String = 'Resume';
        case 'Resume'
            hObject.String = 'Pause';
    end
end
function reset_scale_callback(hObject, eventdata, mechguiopt)

    for i = 1:2
        yyaxis(mechguiopt.ax(i), 'left')
        mechguiopt.ax(i).YLim = [-eps,eps];
        
        yyaxis(mechguiopt.ax(i), 'right')
        mechguiopt.ax(i).YLim = [-eps,eps];
    end
        
    lE = -eps*ones(1,mechguiopt.Numstacks);
    uE = eps*ones(1,mechguiopt.Numstacks);

    yP = [lE, fliplr(uE)];

    for j = 1:numel(mechguiopt.patches)
        mechguiopt.patches(j).YData = yP;
    end
end
function update_frequency_callback(hObject, eventdata, mechguiopt)
    options = hObject.String;

    state = options{hObject.Value};

    if state == "inf"
            mechguiopt.draw_period = Time(0);
    else
        unit = 'kHz';
        num = state(1:end-4);
        num = str2double(num);
        mechguiopt.draw_period = 1 / Frequency(num, unit);
    end
end
function abort_callback(hObject, eventdata, mechguiopt)

    hObject.String = 'Quit';
end
