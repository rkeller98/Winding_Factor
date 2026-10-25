function [fig, plot_tab, method_tab] = investigation_tab(set_id, method, plot_name)
% One reusable figure per set, with method tabs and diagnostic subtabs.
if isnumeric(set_id)
    name = sprintf('Set%d',set_id);
else
    name = char(set_id);
end
tag = ['investigation_' matlab.lang.makeValidName(name)];
fig = findall(groot,'Type','figure','Tag',tag);
if isempty(fig)
    fig = figure('Name',[name ': k-factor investigation'],'Tag',tag, ...
        'Color','w','Visible','on','Position',[30 30 1500 900]);
    theme(fig,'light');
    uitabgroup(fig,'Tag','investigation_methods');
end
tabs = findall(fig,'Type','uitabgroup','Tag','investigation_methods');
method_tag = matlab.lang.makeValidName(method);
method_tab = findall(tabs,'Type','uitab','Tag',method_tag);
if isempty(method_tab)
    method_tab = uitab(tabs,'Title',method,'Tag',method_tag);
    uitabgroup(method_tab,'Tag',[method_tag '_plots']);
end
plots = findall(method_tab,'Type','uitabgroup','Tag',[method_tag '_plots']);
plot_tag = [method_tag '_' matlab.lang.makeValidName(plot_name)];
plot_tab = findall(plots,'Type','uitab','Tag',plot_tag);
if isempty(plot_tab)
    plot_tab = uitab(plots,'Title',plot_name,'Tag',plot_tag);
else
    delete(plot_tab.Children); % Re-running replaces only this diagnostic.
end
tabs.SelectedTab = method_tab;
plots.SelectedTab = plot_tab;
end
