function [mechguiopt] = gui_display(x,signal,time,start_t,t_step,signals,resFiles,mechopt,mnaopt)
%% input data
Numstacks= mechopt.Numstacks;
TWa = rand(size(signal));


%mechfile = 'fBM'
mechfile = 'fMech';

%% GUI init

y_label = ["Amplitude", "displacement"," [nm]","velocity", " [nm/s ?]"];
mechguiopt = mech.gui_v2.init(x, Numstacks, time, signal, NaN, TWa, {'Dropdown-up', 'Dropdown-down'},y_label);
mechguiopt.start_x = 0;
mechguiopt.start = true;
mechguiopt.abort = false;
mechguiopt.subsignals = [];
mechguiopt.global_maximum_x = 0;
mechguiopt.global_minimum_x = 0;
mechguiopt.global_maximum_v = 0;
mechguiopt.global_minimum_v = 0;

fig = mechguiopt.H_FIG;

%% GUI interaction

popup_menu = uicontrol(fig,'Style','popupmenu','Units','normalized');
popup_menu.Position = [0.01 0.85 0.09 0.05];
popup_menu.String = fieldnames(resFiles);

popup_menu_relative = uicontrol(fig,'Style','popupmenu','Units','normalized');
popup_menu_relative.Position = [0.01 0.75 0.09 0.05];
popup_menu_relative.String = ["Absolute","Relative"];

popup_menu_1_1 = uicontrol(fig,'Style','popupmenu','Units','normalized');
popup_menu_1_1.Position = [0.01 0.65 0.09 0.05];
popup_menu_1_1.String = [];

popup_menu_1_2 = uicontrol(fig,'Style','popupmenu','Units','normalized');
popup_menu_1_2.Position = [0.01 0.55 0.09 0.05];
popup_menu_1_2.String = [];

popup_menu_2_1 = uicontrol(fig,'Style','popupmenu','Units','normalized');
popup_menu_2_1.Position = [0.01 0.35 0.09 0.05];
popup_menu_2_1.String = [];

popup_menu_2_2 = uicontrol(fig,'Style','popupmenu','Units','normalized');
popup_menu_2_2.Position = [0.01 0.25 0.09 0.05];
popup_menu_2_2.String = [];

btn_start_pause = uicontrol(fig,'Style','pushbutton','Units', ...
    'normalized','Callback',{@start_stop,popup_menu,mechguiopt});

btn_start_pause.Position = [0 0.92 0.08 0.08];
btn_start_pause.String = {'Start'};

abort = uicontrol(fig,'Style','pushbutton','Units','normalized'...
    ,'Callback',{@abort_fcn,mechguiopt});

abort.Position = [0.90 0 0.08 0.07];
abort.String = {'Abort'};

set(mechguiopt.H_FIG, 'WindowKeyPressFcn', @KeyPress);
mechguiopt.H_FIG.UserData = mechguiopt;


uiwait(fig);

%% Initial conditions

time_value = time.s';
step = round(t_step/(time_value(end)*1e3 / length(time_value)));

[ d, i_start_t ] = min( abs(time_value-start_t));
[row_num,coll_num] = size(resFiles.(mechfile).BMv);
i = i_start_t;

%% Display data
while(~fig.UserData.abort)
    fig.UserData.global_maximum_x = 0;
    fig.UserData.global_minimum_x = 0;
    fig.UserData.global_maximum_v = 0;
    fig.UserData.global_minimum_v = 0;
    while (i < row_num)
        maximum = double.empty(2,0);
        minimum = double.empty(2,0);
        maximum_2 = double.empty(2,0);
        minimum_2 = double.empty(2,0);
        signals_to_update = {};
        
        if(length(fig.UserData.subsignals) > 1)
            if(popup_menu_relative.Value == 1)
                if(isequal(popup_menu.String{popup_menu.Value}, 'fVolt'))
                    [ output_1, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_1_1.Value}, mnaopt, 'transformed','tpoint', i);
                    signals_to_update{1} = output_1;
                    
                    [ output_2, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_1_2.Value}, mnaopt, 'transformed','tpoint', i);
                    signals_to_update{2} = output_2;
                    
                    [ output_3, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_2_1.Value}, mnaopt, 'transformed','tpoint', i);
                    signals_to_update{3} = output_3;
                    
                    [ output_4, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_2_2.Value}, mnaopt, 'transformed','tpoint', i);
                    signals_to_update{4} = output_4;
                else
                    signals_to_update{1} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_1_1.Value})(i,:);
                    signals_to_update{2} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_1_2.Value})(i,:);
                    signals_to_update{3} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_2_1.Value})(i,:);
                    signals_to_update{4} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_2_2.Value})(i,:);
                end
            else
                if(isequal(popup_menu.String{popup_menu.Value}, 'fVolt'))
                    [ output_1, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_1_1.Value}, mnaopt, 'transformed', 'tpoint', i);
                    
                    [ output_1_0, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_1_1.Value}, mnaopt, 'transformed', 'tpoint', 1);
                    
                    signals_to_update{1} = output_1  - output_1_0;
                    
                    [ output_2, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_1_2.Value}, mnaopt, 'transformed','tpoint', i);
                    
                    [ output_2_0, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_1_2.Value}, mnaopt, 'transformed','tpoint', 1);
                    
                    signals_to_update{2} = output_2  - output_2_0;
                    
                    [ output_3, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_2_1.Value}, mnaopt, 'transformed','tpoint', i);
                    
                    [ output_3_0, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_2_1.Value}, mnaopt, 'transformed','tpoint', 1);
                    
                    signals_to_update{3} = output_3  - output_3_0;
                    
                    [ output_4, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_2_2.Value}, mnaopt, 'transformed','tpoint', i);
                    
                    [ output_4_0, ~, ~, ~, ~]  =  getVoltage( resFiles.(popup_menu.String{popup_menu.Value}),...
                        fig.UserData.subsignals{popup_menu_2_2.Value}, mnaopt, 'transformed','tpoint', 1);
                    
                    signals_to_update{4} = output_4 - output_4_0;
                    
                    if i == 0
                        
                    end
                    
                else
                    signals_to_update{1} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_1_1.Value})(i,:)...
                        - resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_1_1.Value})(1,:);
                    
                    signals_to_update{2} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_1_2.Value})(i,:)...
                        - resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_1_2.Value})(1,:);
                    
                    signals_to_update{3} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_2_1.Value})(i,:)...
                        - resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_2_1.Value})(1,:);
                    
                    signals_to_update{4} = resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_2_2.Value})(i,:)...
                        - resFiles.(popup_menu.String{popup_menu.Value}).(fig.UserData.subsignals{popup_menu_2_2.Value})(1,:);
                end
            end
            
            % BMx
            [maximum(1),minimum(1)] = update_line( ...
                mechguiopt.displacement(1), signals_to_update{1});
            % TMx
            [maximum(2),minimum(2)] = update_line( ...
                mechguiopt.displacement(2),signals_to_update{2});
            
            if(max(maximum) > fig.UserData.global_maximum_x)
                fig.UserData.global_maximum_x = max(maximum);
            end
            
            if(min(minimum) < fig.UserData.global_minimum_x)
                fig.UserData.global_minimum_x = min(minimum);
            end
            
            update_ylim(mechguiopt.ax_displacement, fig.UserData.global_maximum_x,fig.UserData.global_minimum_x);
            
            % BMv
            [maximum_2(1),minimum_2(1)] = update_line( ...
                mechguiopt.velocity(1),signals_to_update{3});
            % TMv
            [maximum_2(2),minimum_2(2)] = update_line( ...
                mechguiopt.velocity(2), signals_to_update{4});
            
            if(max(maximum_2) > fig.UserData.global_maximum_v)
                fig.UserData.global_maximum_v = max(maximum_2);
            end
            
            if(min(minimum_2) < fig.UserData.global_minimum_v)
                fig.UserData.global_minimum_v = min(minimum_2);
            end
            
            update_ylim(mechguiopt.ax_velocity,fig.UserData.global_maximum_v,fig.UserData.global_minimum_v)
            
            set(mechguiopt.line_0, 'Xdata', [time_value(i)*1e3, time_value(i)*1e3], ...'Ydata', [-maxsig, maxsig],...
                'Color', 'r');
            drawnow;
            
            if(mechguiopt.H_FIG.UserData.start_x > 0)
                i = round(mechguiopt.H_FIG.UserData.start_x/(time_value(end)*1e3 / length(time_value)));
                mechguiopt.H_FIG.UserData.start_x = 0;
                continue
            end
            i = i + step;
            if(fig.UserData.abort)
                return
            end
        else
            uiwait(fig);
            if(mechguiopt.H_FIG.UserData.start_x > 0)
                i = round(mechguiopt.H_FIG.UserData.start_x/(time_value(end)*1e3 / length(time_value)));
                mechguiopt.H_FIG.UserData.start_x = 0;
                continue
            end
        end
        
    end
    i = i_start_t;
    fig.UserData.ax_displacement.YLim = 1.1*[min(min(signals_to_update{1}),min(signals_to_update{2})),max(max(signals_to_update{1}),max(signals_to_update{2}))];
    fig.UserData.ax_velocity.YLim = 1.1*[min(min(signals_to_update{3}),min(signals_to_update{4})),max(max(signals_to_update{3}),max(signals_to_update{4}))];
end

%%
    function [maximum,minimum] = update_line(line_handle, Y)
        line_handle.YData = Y;
        maximum = max(Y);
        minimum = min(Y);
    end

    function update_ylim(ax_handle, maximum,minimum)
        if(max(maximum) == 0)
            maximum = 1.0e-20*0.2220;
        end
        if(min(minimum) == 0)
            minimum = 1.0e-20*-0.2220;
        end
        %ax_handle.YLim =  [min(minimum) - 0.1*(max(maximum)-min(minimum)),max(maximum) + 0.1*(max(maximum)-min(minimum))];
        ax_handle.YLim =  [min(minimum),max(maximum)];
    end

    function KeyPress(Source, EventData)
        if(EventData.Character == 't')
            [time_x, y] = selectDataPoints(Source.UserData.ax_sig);
            Source.UserData.start_x = time_x;
            Source.UserData.ax_displacement.YLim = 1.0e-20*[-0.2220,0.2220];
            Source.UserData.ax_velocity.YLim = 1.0e-20*[-0.2220,0.2220];
        end
        
    end

    function [x, y] = selectDataPoints(ax)
        roi = drawpoint(ax);
        x = roi.Position(1);
        y = roi.Position(2);
    end

    function start_stop(btn,event,popup,source)
        source.H_FIG.UserData.start = ~source.H_FIG.UserData.start;
        sub_signals = fieldnames(resFiles.(popup.String{popup.Value}));
        sub_signal = 1;
        source.H_FIG.UserData.global_maximum_x = 0;
        source.H_FIG.UserData.global_minimum_x = 0;
        source.H_FIG.UserData.global_maximum_v = 0;
        source.H_FIG.UserData.global_minimum_v = 0;
        
        
        
        while sub_signal < numel(sub_signals)
            [rows, colls] = size(resFiles.(popup.String{popup.Value}).(sub_signals{sub_signal}));
            if(colls ~= numel(x) || rows == 1)
                sub_signals(sub_signal) = [];
                continue
            end
            sub_signal = sub_signal + 1;
        end
        if(isequal(popup_menu.String{popup_menu.Value}, 'fVolt'))
            names_to_add = {'ScalaMedia','IHC_ic','OHC_ic','IHC_ec','OHC_ec','SpiralLigament','IHC','OHC'}';
            sub_signals = [names_to_add;sub_signals];
        elseif(isequal(popup_menu.String{popup_menu.Value}, 'fCurr'))
            names_to_add = {'IHC','OHC','IHC_MET','OHC_MET','StV'}';
            sub_signals = [names_to_add;sub_signals];
        end
        source.H_FIG.UserData.subsignals = sub_signals;
        popup_menu_1_1.String = sub_signals;
        popup_menu_1_2.String = sub_signals;
        popup_menu_2_1.String = sub_signals;
        popup_menu_2_2.String = sub_signals;
        
        if (popup.Value == 1)
            source.ax_displacement.YLabel.String = 'displacement [nm]';
            source.ax_velocity.YLabel.String = 'velocity [nm/s ?]';
        else
            source.ax_displacement.YLabel.String = 'Test';
            source.ax_velocity.YLabel.String = 'Test';
        end
        
        if(source.H_FIG.UserData.start)
            btn.String = {'Start'};
            uiwait(fig);
            
        else
            btn.String = {'Pause'};
            uiresume(fig);
        end
    end

    function abort_fcn(btn,event,source)
        names = fieldnames(source);
        source.H_FIG.UserData.abort = ~source.H_FIG.UserData.abort;
    end
end


