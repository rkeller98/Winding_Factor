function export_investigation_plot(target, filename, export_format)
% Export only the requested axes/layout, activating its tabs temporarily.
if isempty(export_format)
    return;
end
containers = gobjects(0);
selected = gobjects(0);
item = target;
while ~isa(item,'matlab.ui.Figure')
    if isa(item,'matlab.ui.container.Tab')
        containers(end+1) = item.Parent;
        selected(end+1) = item.Parent.SelectedTab;
        item.Parent.SelectedTab = item;
    end
    item = item.Parent;
end
drawnow;
export_analysis_plot(target,filename,export_format,150);
for j = 1:numel(containers)
    containers(j).SelectedTab = selected(j);
end
end
