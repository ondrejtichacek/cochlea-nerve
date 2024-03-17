function status = update( t, y, flag, mechguiopt, mechopt, mnaopt )
%UPDATE

mechguiopt.paused = false;

if isempty(flag) || strcmp(flag, 'loader')
    while strcmp(mechguiopt.uicontrol.btn_start.String, 'Resume')
        mechguiopt.paused = true;
        pause(0.1)
        if strcmp(flag, 'loader')
            newflag = 'loader';
        else
            newflag = 'update';
        end
        % mech.gui_v2.update(t, y, newflag, mechguiopt, mechopt, mnaopt);
        status = false;
        return
    end
    if strcmp(mechguiopt.uicontrol.btn_abort.String, 'Quit')
        mechguiopt.aborted = true;
        status = true;
        return
    end
    
    if abs(t(end) - mechguiopt.last_draw) < mechguiopt.draw_period.s
        status = false;
        return
    end
    
end

if isempty(flag) || strcmp(flag, 'update')
    t = t(end);
    y = y(:,end);
end
if strcmp(flag, 'loader')
    t_index = mechguiopt.loader_t_index;
else
    t_index = 1;
end

channels = {mnaopt.channels.name};

channel_signals_p = {};
channel_signals_dp = {};

for i = 1:numel(channels)
    channel_signals_p{end+1} = ['p_channel_', channels{i}];
    channel_signals_dp{end+1} = ['dp_channel_', channels{i}];
end

%%

if strcmp(flag, 'loader') && mechguiopt.preloaded == false
    
    dbprintf('Preloading synapse results ...\n')

    variables = {'release_event', 'action_potential'};    
    
    for i = 1:numel(variables)
        variable = variables{i};

        if ~isempty(y.synapse) || ~isempty(y.nerve)
            switch variable
                case 'release_event'
                    TT = {};
                    XX = {};
                    for l = 1:numel(y.synapse)
                        tt = y.synapse(l).t.s;
                        R = y.synapse(l).vesicle_release_events;
                        for kkk = 1:numel(R)
                            ttt = R{kkk}' * 1e3;
                            TT{l,kkk} = ttt;
                            % XX{l,kkk} = ones(size(ttt)) * y.synapse(l).xgrid(kkk);
                            XX{l,kkk} = ones(size(ttt)) * y.synapse(l).xgrid_ind(kkk) / mnaopt.Numstacks;
                        end
                    end
                case 'action_potential'
                    TT = {};
                    XX = {};
                    for l = 1:numel(y.nerve)
                        tt = y.nerve(l).t.s;
                        FS = 1/mean(diff(tt));
                        VV = y.nerve(l).V;
                        for kkk = 1:size(VV, 2)
                            [pks, ttt] = findpeaks(VV(:,kkk), FS, ...
                                'MinPeakHeight', 100, ...
                                'MinPeakDistance', 0.001);
                            TT{l,kkk} = ttt * 1e3;
                            % XX{l,kkk} = ones(size(ttt)) * y.nerve(l).xgrid(kkk);
                            XX{l,kkk} = ones(size(ttt)) * y.nerve(l).xgrid_ind(kkk) / mnaopt.Numstacks;
                        end
                    end
                otherwise
                    error('not set up')
            end
            
            bs = Time(0.5, 'ms');
            nbins = round(tt(end) / bs.s);
            counts = zeros(size(y.nerve(l).V, 2), nbins);
            for kkk = 1:size(y.nerve(l).V, 2)
                [counts(kkk,:), edges] = histcounts(cat(1,TT{:,kkk}), nbins);
            end

            mechguiopt.preloaded_vars.(variable).events = TT;
            mechguiopt.preloaded_vars.(variable).xpos = XX;
            mechguiopt.preloaded_vars.(variable).counts = counts;
            mechguiopt.preloaded_vars.(variable).edges = edges;
        end

    end
    dbprintf('Preloading synapse results done\n')

    dbprintf('Preloading other results ...\n')
    
    name_part = 'IHC';
           
    val = getVoltage(y.oc_mna.fVolt, name_part, mnaopt, 'transformed');
    tt = y.oc_mna.fT.T;
    
    mechguiopt.preloaded_vars.(sprintf('V_%s', name_part)).V = val;
    mechguiopt.preloaded_vars.(sprintf('V_%s', name_part)).t = tt;

    dbprintf('Preloading other results done\n')
    mechguiopt.preloaded = true;
end
%%

if isempty(flag) || strcmp(flag, 'update') || strcmp(flag, 'loader')

    % dbprintf('h = %s\n', Time(t - mechguiopt.last_t, 's'));
    mechguiopt.last_t = t;
    
    for kk = 1:2
        for ll = 1:2
            for k = 1:numel(mechguiopt.uicontrol.signal(kk,:,ll))
                
                S = mechguiopt.uicontrol.signal(kk,k,ll);
            
                val = S.Value;
                str = S.String;
                sel = str{val};

                mechguiopt.text_label(kk,k,ll).String = get_display_name(sel);
                
                sig = get_signal(sel, y, t_index);
                
                R = mechguiopt.uicontrol.rel_switch(kk,k,ll);
            
                val = R.Value;
                str = R.String;
                rel = str{val};
                
                if strcmp(rel, 'Relative')
                    if strcmp(flag, 'loader')
                        y0 = y;
                        t_index_0 = 1;
                    else
                        y0 = mechguiopt.y0;
                        t_index_0 = 1;
                    end
                    rel_sig = get_signal(sel, y0, t_index_0);
                    sig = sig - rel_sig;
                end
                
                m(kk,k,ll) = update_line( ...
                        mechguiopt.lines(kk,k,ll), ...
                        sig, sel);
                
                if mechguiopt.draw_patches
                    update_patch( ...
                            mechguiopt.patches(kk,k,ll), ...
                            sig, sel)
                end
            end
            if ll == 1
                yyaxis(mechguiopt.ax(kk), 'left')
            elseif ll == 2
                yyaxis(mechguiopt.ax(kk), 'right')
            end
            update_ylim(mechguiopt.ax(kk), max(m(kk,:,ll)));
        end
    end

    for kk = 3
        % axes(mechguiopt.ax(kk))

        % V = y.fMech.BMx';
        % V = V(:,ind);

        name_part = 'IHC';
            
        if strcmp(flag, 'loader')
            % ind = 1:10:t_index;
            % val = getVoltage(y.oc_mna.fVolt, name_part, mnaopt, 'transformed', ind);
            % tt = y.oc_mna.fT.T;
            % tt = tt(ind);

            ind = 1:t_index;
            val = mechguiopt.preloaded_vars.(sprintf('V_%s', name_part)).V;
            tt = mechguiopt.preloaded_vars.(sprintf('V_%s', name_part)).t;

            tt = tt(ind);
            val = val(ind,:);
        else
            val = getVoltage(y', name_part, mnaopt, 'original');

            if isempty(mechguiopt.ax(kk).UserData)
                mechguiopt.ax(kk).UserData.tt = [];
                mechguiopt.ax(kk).UserData.val = [];
            end
            
            tt = [mechguiopt.ax(kk).UserData.tt, t];
            val = [mechguiopt.ax(kk).UserData.val; val];

            mechguiopt.ax(kk).UserData.tt = tt;
            mechguiopt.ax(kk).UserData.val = val;
            
        end

        xx = tt*1e3;

        V = val';
        V = V - V(:,1);

        % yy = mnaopt.xgrid;
        % imagesc('XData', xx, 'YData', yy, 'CData', V, [0,1e-6])

        % M = max(abs(V(:)));
        M = max(V(:));

        % mechguiopt.ax(kk).CLim = [-M-1e3*eps, M+1e3*eps];
        mechguiopt.ax(kk).CLim = [0, M+1e3*eps];

        % [XX, YY] = meshgrid(xx, mnaopt.xgrid);

        set(mechguiopt.ax(kk).Children(4), ... 
            'CData', V, ...
            ...'XData', XX, ...
            ...'YData', YY)
            'XData', xx, ...
            'YData', mnaopt.xgrid);

        % xlim([tt(1), tt(end)])
        
        mock_action_potentials = false;
        if mock_action_potentials
            fac = mechguiopt.draw_period / Time(0.01, 'ms');    
            scale = linf(3e-2, 1e-2, mnaopt.xgrid);
            P = V(:,end) ./ scale(:) * fac;
    
            xxx = rand(mnaopt.Numstacks,1) < P;
            if any(xxx)
                nt = sum(xxx);
                ttt = mechguiopt.ax(kk).Children(2).XData;
                ttt = [ttt, repmat(t*1e3,1,nt)];
                yyy = mechguiopt.ax(kk).Children(2).YData;
                
                % do it like this because imagesc is always uniform
                % y_new = mnaopt.xgrid(xxx);
                y_new = linspace(0,1,mnaopt.Numstacks);
                y_new = y_new(xxx);
                
                yyy = [yyy, y_new];
                mechguiopt.ax(kk).Children(2).XData = ttt;
                mechguiopt.ax(kk).Children(2).YData = yyy;
            end
        end

        S = mechguiopt.uicontrol.heat_over_drop(1,1,1);
        val = S.Value;
        str = S.String;
        sel = str{val};

        % variable = 'release_event';
        % variable = 'action_potential';
        variable = sel;

        if strcmp(variable, 'select')
            set(mechguiopt.ax(kk).Children(2), ...
                'XData', [], ...
                'YData', []);
        else
            if ~isemptystruct(mechguiopt.preloaded_vars)
                X = cat(1,mechguiopt.preloaded_vars.(variable).events{:});
                Y = cat(1,mechguiopt.preloaded_vars.(variable).xpos{:});
                ind = X <= t*1e3;
                set(mechguiopt.ax(kk).Children(2), ...
                    'XData', X(ind), ...
                    'YData', Y(ind));
            end
        end

    end
    
    for kk = 4
        if ~isemptystruct(mechguiopt.preloaded_vars)
            variable = 'release_event';
            variable = 'action_potential';
    
            counts = mechguiopt.preloaded_vars.(variable).counts;
            C = counts .* (counts > 0.1*numel(y.nerve));
            i_activity = C > 0;
            C(i_activity) = 0.5 + C(i_activity)/2;
            A = (counts > 0.1*numel(y.nerve));

            edges = mechguiopt.preloaded_vars.(variable).edges;
            centers = (edges(1:end-1) + edges(2:end))/2;
            t_ind = find(t*1e3 <= centers, 1);
            if isempty(t_ind)
                t_ind = 1;
            end
    
            set(mechguiopt.ax(kk).Children(3), ... 
                'CData', C(:,1:t_ind), ...
                ...'AlphaData', A(:,1:t_ind), ...
                ...'XData', XX, ...
                ...'YData', YY)
                'XData', centers(1:t_ind))%, ...
                ...'YData', mnaopt.xgrid);
            
    
            S = mechguiopt.uicontrol.heat_over_drop(2,1,1);
            val = S.Value;
            str = S.String;
            sel = str{val};
    
            % variable = 'release_event';
            % variable = 'action_potential';
            variable = sel;
        end

        if strcmp(variable, 'select')
            set(mechguiopt.ax(kk).Children(2), ...
                'XData', [], ...
                'YData', []);
        else
            if ~isemptystruct(mechguiopt.preloaded_vars)
                X = cat(1,mechguiopt.preloaded_vars.(variable).events{:});
                Y = cat(1,mechguiopt.preloaded_vars.(variable).xpos{:});
                ind = X <= t*1e3;
                set(mechguiopt.ax(kk).Children(2), ...
                    'XData', X(ind), ...
                    'YData', Y(ind));
            end
        end
    end

    % time
    set(mechguiopt.line_0, 'Xdata', [t*1e3, t*1e3], ...'Ydata', [-maxsig, maxsig],...
        'Color', 'r');

    drawnow;
    mechguiopt.last_draw = t;
    
elseif strcmp(flag, 'init')
    % init
    
    t = t(end);
    y = y(:,end);
    
    mechguiopt.y0 = y;
    
elseif strcmp(flag, 'done')
    % cleanup
else
    error('Unrecognized flag value `%s`', disp(flag));
end

status = false;

    function m = update_line(line_handle, Y, name)
        
        line_handle.YData = Y;
        % line_handle.DisplayName = name;
        
        m = max(abs(Y));
    end
    function update_patch(patch_handle, Y, name)
        if name == "select"
            patch_handle.YData = [];
        else
            if isempty(patch_handle.YData)
                patch_handle.YData = [inf(mechguiopt.Numstacks,1); -inf(mechguiopt.Numstacks,1)];
            end
            patch_handle.YData = [ ...
                min(patch_handle.YData(1:mechguiopt.Numstacks), Y);
                max(patch_handle.YData(mechguiopt.Numstacks+1:end), flip(Y))];
        end
    end
    function update_ylim(ax_handle, M)
        if M > ax_handle.YLim(2)
            % disp(ax_handle.YAxisLocation)
            ax_handle.YLim = 1.05*[-M, M];
        end
    end
    function val = get_display_name(name)
        name = char(name);

        display_names = struct(...
            "mech_BMx", "BM displacement", ...
            "mech_BMv", "BM velocity", ...
            "mech_TMx", "TM displacement", ...
            "mech_TMv", "TM velocity", ...
            "volt_IHC", "IHC receptor potential", ...
            "volt_OHC", "OHC receptor potential", ...
            "select", "");

        if isfield(display_names, name)
            val = display_names.(name);
        else
            warning('No display name configured for signal %s', name);
            val = '';
        end
    end
    function val = get_signal(name, y, t_index)
        name = char(name);

        if name == "select"
            val = NaN(mechguiopt.Numstacks, 1);
                
        elseif any(name == ["mech_BMx", "mech_BMv", "mech_TMx", "mech_TMv"])
            
            if strcmp(flag, 'loader')
                p = name(6:8);
                val = y.oc_mna.fMech.(p)(t_index,:);
            else
                p = name(6:7);
                if p == "BM"
                    K = mechopt.NMkonst_BM;
                elseif p == "TM"
                    K = mechopt.NMkonst_OHC_cilia;
                end
                
                i = mechguiopt.size_info.(name).start;
                j = mechguiopt.size_info.(name).end;
    
                val = K*y(i:j,end);
            end
                
        
        elseif any(name == ["stapes_mech_BMx", "stapes_mech_TMx", "stapes_mech_BMv", "stapes_mech_TMv"])
            
            if strcmp(flag, 'loader')
                warning('Signal %s not accessible in ''loader'' mode', name)
                val = NaN(mechguiopt.Numstacks, 1);
            else
                np = name(numel('stapes_')+1:end);
                
                i = mechguiopt.size_info.mech.start;
                j = mechguiopt.size_info.mech.end;
                    
                g = mechopt.OdeFun_t;
                val = g(t);
                                    
                ii = mechguiopt.size_info.(np).start;
                jj = mechguiopt.size_info.(np).end;
    
                ii = ii - i + 1;
                jj = jj - i + 1;
    
                val = val(ii:jj);
            end
        
        elseif any(name == ["amplifier_mech_BMv", "amplifier_mech_TMv"])
            
            if strcmp(flag, 'loader')
                warning('Signal %s not accessible in ''loader'' mode', name)
                val = NaN(mechguiopt.Numstacks, 1);
            else
                np = name(numel('amplifier_')+1:end);
                
                i = mechguiopt.size_info.mech.start;
                j = mechguiopt.size_info.mech.end;

                if mechopt.amplifier == "electric"
                    vohc = get_signal("volt_OHC", y, t_index);
                    h = mechopt.OdeFun_v;
                    val = h(vohc - mechopt.vohc_0);

                elseif mechopt.amplifier == "mechanic"
        
                    h = mechopt.OdeFun_y;
                    val = h(y(i:j,end));
                end
                    
                ii = mechguiopt.size_info.(np).start;
                jj = mechguiopt.size_info.(np).end;
    
                ii = ii - i + 1;
                jj = jj - i + 1;
    
                val = val(ii:jj);
            end

        elseif any(name == ["F_OHC"])
            
            if strcmp(flag, 'loader')
                error('not yet implemented')
            else
                np = 'dvohc';
                    
                i = mechguiopt.size_info.(np).start;
                j = mechguiopt.size_info.(np).end;
    
                val = y(i:j);
            end
                
        elseif name == "amplification_BM"
            BMv = get_signal("mech_BMv", y, t_index);
            ampl = get_signal("amplifier_mech_BMv", y, t_index);

            val = ampl .* BMv ;
            % val = sign(ampl .* BMv) .* abs(ampl) / max(abs(BMv));
                
        elseif name == "amplification_TM"
            TMv = get_signal("mech_TMv", y, t_index);
            ampl = get_signal("amplifier_mech_TMv", y, t_index);

            val = TMv .* ampl;
            
        elseif any(name == ["curr_IHC", "curr_OHC"])
         
            name_part = name(numel('curr_')+1:end);

            val = getCurrent(y', name_part, mnaopt);
        
        elseif any(name == ["tol_curr_IHC", "tol_curr_OHC"])
            
            name_part = name(numel('tol_curr_')+1:end);

            [~, v_sel] = getCurrent(y', name_part, mnaopt);

            ii = mechguiopt.size_info.circuit.start;
            jj = mechguiopt.size_info.circuit.end;

            if isa(mechguiopt.abs_tol, 'function_handle')
                val = mechguiopt.abs_tol(t);
                val = val(ii:jj);
            else
                val = mechguiopt.abs_tol(ii:jj);
            end

            val = val(v_sel);

        elseif any(name == ["volt_IHC", "volt_OHC", ...
                            "volt_IHC_ic", "volt_OHC_ic", ...
                            "volt_IHC_ec", "volt_OHC_ec"])
         
            name_part = name(numel('volt_')+1:end);
            
            if strcmp(flag, 'loader')
                val = getVoltage(y.oc_mna.fVolt, name_part, mnaopt, 'transformed', t_index);
            else
                val = getVoltage(y', name_part, mnaopt, 'original');
            end
        
        elseif any(strcmp(name, channel_signals_p))
            
            name_part = name(numel('p_')+1:end);

            ii = mechguiopt.size_info.(name_part).start;
            jj = mechguiopt.size_info.(name_part).end;

            val = y(ii:2:jj);

        elseif any(strcmp(name, channel_signals_dp))
            
            name_part = name(numel('dp_')+1:end);

            ii = mechguiopt.size_info.(name_part).start;
            jj = mechguiopt.size_info.(name_part).end;

            ii = ii + 1;

            val = y(ii:2:jj);

        elseif any(name == ["tol_volt_IHC_ic", "tol_volt_IHC_ec", "tol_volt_OHC_ic", "tol_volt_OHC_ec"])
            
            name_part = name(numel('tol_volt_')+1:end);

            [~, v_sel, ~, ~, ~] = getVoltage(y', name_part, mnaopt, 'original');

            ii = mechguiopt.size_info.circuit.start;
            jj = mechguiopt.size_info.circuit.end;

            if isa(mechguiopt.abs_tol, 'function_handle')
                val = mechguiopt.abs_tol(t);
                val = val(ii:jj);
            else
                val = mechguiopt.abs_tol(ii:jj);
            end

            val = val(v_sel);

        elseif any(name == ["tol_mech_BMx", "tol_mech_TMx", "tol_mech_BMv", "tol_mech_TMv"])
        
            name_part = name(numel('tol_')+1:end);
            p = name(numel('tol_mech_')+1:numel('tol_mech_')+2);
            
            if p == "BM"
                K = mechopt.NMkonst_BM;
            elseif p == "TM"
                K = mechopt.NMkonst_OHC_cilia;
            end
            
            ii = mechguiopt.size_info.(name_part).start;
            jj = mechguiopt.size_info.(name_part).end;

            if isa(mechguiopt.abs_tol, 'function_handle')
                val = mechguiopt.abs_tol(t);
                val = val(ii:jj);
            else
                val = mechguiopt.abs_tol(ii:jj);
            end
            val = K * val(:);

        elseif regexp(name, 'R_*')
        
            name_part = name(numel('R_')+1:end);
            
            % R = mnaopt.circuit.get_element_by_name(name_part);

            if strcmp(flag, 'loader')
                val = getVoltage(y.oc_mna.fVolt, name_part, mnaopt, 'transformed', t_index);
            else
                val = getVoltage(y', name_part, mnaopt, 'original');
            end
            
        else
            error('unknown option %s', name)
        end

        val = val(:);

    end

end
