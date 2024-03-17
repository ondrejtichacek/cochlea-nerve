function [] = loader(results, stimulus, mechopt, mnaopt)
%LOADER 

fT = results.oc_mna.fT;

t = fT.T;

t0 = Time(0, 'ms');

mechguiopt = mech.gui_v2.init( ...
    mnaopt.xgrid, stimulus.audiotime, stimulus.audio, mnaopt, ...
    t0=t0, ...
    init_state='paused');

mechguiopt.draw_period = Time(0.01, 'ms');
% mechguiopt.draw_period = Time(1, 'ms');
% mechguiopt.y0 = resFiles;

% preload_res_files = false;
preload_res_files = true;

% parts = fieldnames(results);
dbprintf('Preloading results ...\n')
% parts = {'oc_mna', 'middle_ear'}; % preloading of other parts not supported
parts = {'oc_mna'}; % preloading of other parts not supported
if preload_res_files
    for j = 1:numel(parts)
        resFiles = results.(parts{j});

        num_replicas = numel(resFiles);
        fn = fieldnames(resFiles);

        loadedFiles = struct();
        
        for k = num_replicas:-1:1
            for i = 1:numel(fn)
                f = resFiles(k).(fn{i}).Properties.Source;
                if exist(f, "file")
                    loadedFiles(k).(fn{i}) = load(f);
                end
            end
        end
        results.(parts{j}) = loadedFiles;
    end
end
dbprintf('Preloading done.\n')

i = 1;
i = find(t0.s <= t, 1);

while true
    
    mechguiopt.loader_t_index = i;
    
    do_break = mech.gui_v2.update(t(i), results, 'loader', mechguiopt, mechopt, mnaopt);
    if do_break == true
        break
    end

    if ~isempty(mechguiopt.loader_t_click)
        
        t_click = mechguiopt.loader_t_click;
        mechguiopt.loader_t_click = [];
        
        ind = find(t > t_click, 1);
        if ~isempty(ind)
            i = ind - 1;
        else
            warining('Clik out of range')
        end
    end

    if ~mechguiopt.paused
        i = i + 1;
        if i > numel(t)
            break
        end
    end
end

% for i = 1:numel(t)
%     mechguiopt.loader_t_index = i;
%     do_break = mech.gui_v2.update(t(i), results, 'loader', mechguiopt, mechopt, mnaopt);
%     if do_break == true
%         break
%     end
% end
end

