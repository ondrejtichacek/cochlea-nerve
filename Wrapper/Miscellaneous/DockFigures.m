function [dockedfigures, hgroup] = DockFigures(dock_group_spec, dockedfigures, hgroup)
%DOCKFIGURES
arguments
    dock_group_spec
    dockedfigures = []
    hgroup = struct()
end

try
    ff = findall(groot,'Type','figure');
    for i = 1:numel(ff)
        if ~any(find(dockedfigures == ff(i).Number))
            dockedfigures(end+1) = ff(i).Number;
            hgroup.(dock_group_spec) = setFigDockGroup(ff(i), dock_group_spec);
        end
    end
catch
end

end

