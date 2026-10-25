function fig = plot_local_factor(local, values, field_name, valid, export_format)
% Physical current meshes, separately connected within each nominal group.
% plot_local_factor(local, 'psi_d'/'psi_q'/'ls'/'kd'/'kq')
% plot_local_factor(local, local.kd, 'k_d', local.valid_pair, 'png')
requested_valid = [];
if nargin >= 4
    requested_valid = valid;
end
if nargin < 2 || isempty(values)
    values = local.klocal;
end
if nargin < 3 || isempty(field_name)
    field_name = 'k_{loc}';
end
if nargin < 4 || isempty(valid)
    valid = local.valid_local;
end
if nargin < 5
    export_format = '';
end
plot_name = field_name;
if strcmp(field_name,'k_{loc}')
    plot_name = 'k LS';
end
if ischar(values) || isstring(values)
    switch char(values)
        case 'psi_d'
            values = local.kd_from_psi_d;
            field_name = 'k_d^{(\psi_d)}';
            plot_name = 'k aus psi_d';
            valid = local.valid_psi_d;
        case 'psi_q'
            values = local.kd_from_psi_q;
            field_name = 'k_d^{(\psi_q)}';
            plot_name = 'k aus psi_q';
            valid = local.valid_psi_q;
        case 'ls'
            values = local.klocal;
            field_name = 'k_{LS}';
            plot_name = 'k LS';
            valid = local.valid_local;
        case 'kd'
            values = local.kd;
            field_name = 'k_d (two-axis)';
            plot_name = 'k_d (2-achsig)';
            valid = local.valid_pair;
        case 'kq'
            values = local.kq;
            field_name = 'k_q (two-axis)';
            plot_name = 'k_q (2-achsig)';
            valid = local.valid_pair;
    end
end
if ~isempty(requested_valid)
    valid = valid & requested_valid;
end
x = local.x;
values = values(:);
valid = valid(:) & isfinite(values) & all(isfinite(x), 2);
ng = numel(local.excitation);
scale = max(abs(x(:, 1:2)), [], 1);
scale(scale == 0) = 1;
[fig,plot_tab,method_tab] = investigation_tab(local.set,'Lokale Faktoren',plot_name);
ax = axes(plot_tab, 'Position', [.04 .16 .33 .69]);
heat_ax = axes(plot_tab,'Position',[.48 .16 .25 .69]);
hold(ax, 'on');
meshes = gobjects(ng, 1);
points = gobjects(ng, 1);
checks = gobjects(ng, 1);
panel = uipanel(plot_tab, 'Title', 'Erregerstromgruppen', 'Units', 'normalized', ...
    'Position', [.83 .04 .16 .92]);
uicontrol(panel, 'Style', 'pushbutton', 'String', 'Alle', 'Tag', 'show_all_groups', ...
    'Units', 'normalized', 'Position', [.04 .92 .43 .05], 'Callback', @(~, ~) show_all(true));
uicontrol(panel, 'Style', 'pushbutton', 'String', 'Keine', 'Tag', 'hide_all_groups', ...
    'Units', 'normalized', 'Position', [.53 .92 .43 .05], 'Callback', @(~, ~) show_all(false));
height = min(.05, .65/ng);

for g = 1:ng
    previous = findall(method_tab,'Tag',sprintf('local_group_%02d',g));
    visible = true;
    if ~isempty(previous)
        visible = logical(previous(1).Value);
    end
    rows = find(local.groups == g);
    xy = x(rows, 1:2)./scale;
    % Triangulate ALL original points first. Removing invalid vertices later
    % leaves holes rather than reconnecting the remaining vertices across them.
    dt = delaunayTriangulation(xy);
    tri = dt.ConnectivityList;
    edges = dt.edges;
    limit = 9*mean(sum((xy(edges(:, 1), :)-xy(edges(:, 2), :)).^2, 2));
    keep = all(valid(rows(tri)), 2);
    for edge = [1 2; 1 3; 2 3]'
        keep = keep & sum((xy(tri(:, edge(1)), :)-xy(tri(:, edge(2)), :)).^2, 2) < limit;
    end
    z = values(rows);
    z(~valid(rows)) = nan;
    meshes(g) = patch(ax, 'Faces', tri(keep, :), 'Vertices', [x(rows, 1:2) z], ...
        'FaceVertexCData', x(rows, 3), 'FaceColor', 'none', 'EdgeColor', 'interp', ...
        'LineWidth', .65, 'Tag', sprintf('local_mesh_%02d', g));
    rows = rows(valid(rows));
    points(g) = scatter3(ax, x(rows, 1), x(rows, 2), values(rows), 9, x(rows, 3), ...
        'filled', 'Tag', sprintf('local_points_%02d', g));
    checks(g) = uicontrol(panel, 'Style', 'checkbox', 'String', ...
        sprintf('%g A (Soll)', local.excitation(g)), 'Value', double(visible), 'Units', 'normalized', ...
        'Position', [.04 .87-g*height .94 height], 'Tag', sprintf('local_group_%02d', g), ...
        'Callback', @(source, ~) show_group(g, logical(source.Value)));
    show_group(g,visible);
end
uicontrol(panel,'Style','text','String','Heatmap-Gruppe','Units','normalized', ...
    'Position',[.04 .14 .92 .035]);
uicontrol(panel,'Style','popupmenu','String',compose('%g A (Soll)',local.excitation), ...
    'Units','normalized','Position',[.04 .09 .92 .045],'Value',1, ...
    'Callback',@(source,~) draw_heatmap(source.Value));
draw_heatmap(1);
uicontrol(panel,'Style','text','String',{'Auswahl gilt in allen','Faktor-Tabs.'}, ...
    'Units','normalized','Position',[.04 .01 .92 .07]);
xlabel(ax, 'i_d in A');
ylabel(ax, 'i_q in A');
zlabel(ax, sprintf('%s in A/A', field_name));
title(ax, sprintf('Set%d: %s; saved 80-neighbour derivatives', local.set, field_name));
grid(ax, 'on');
view(ax, 35, 25);
colormap(ax, turbo(256));
cb = colorbar(ax);
cb.Position = [.39 .16 .012 .66];
cb.Label.String = 'Measured i_e in A';
clim(ax, [min(x(:, 3)) max(x(:, 3))]);
rotate3d(fig, 'on');
if ~isempty(export_format)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'matlab', 'common'));
    name = regexprep(field_name, '[^a-zA-Z0-9]', '');
    export_investigation_plot(ax, fullfile(root, 'matlab', 'k_factor_investigation', 'results', ...
        sprintf('local_%s_set%d', name, local.set)), export_format);
    export_investigation_plot(heat_ax, fullfile(root, 'matlab', 'k_factor_investigation', 'results', ...
        sprintf('local_%s_heatmap_set%d', name, local.set)), export_format);
end

    function show_group(g, visible)
        state = 'off';
        if visible
            state = 'on';
        end
        set(findall(method_tab,'Tag',sprintf('local_mesh_%02d',g)),'Visible',state);
        set(findall(method_tab,'Tag',sprintf('local_points_%02d',g)),'Visible',state);
        set(findall(method_tab,'Tag',sprintf('local_heatmap_%02d',g)),'Visible',state);
        set(findall(method_tab,'Tag',sprintf('local_group_%02d',g)),'Value',double(visible));
    end

    function show_all(visible)
        set(checks, 'Value', double(visible));
        for g = 1:ng
            show_group(g, visible);
        end
    end

    function draw_heatmap(g)
        cla(heat_ax);
        rows = local.groups==g & valid;
        h = scatter(heat_ax,x(rows,1),x(rows,2),18,values(rows),'filled', ...
            'Tag',sprintf('local_heatmap_%02d',g));
        check = findall(plot_tab,'Tag',sprintf('local_group_%02d',g));
        if ~check.Value
            h.Visible = 'off';
        end
        xlabel(heat_ax,'i_d in A');
        ylabel(heat_ax,'i_q in A');
        title(heat_ax,sprintf('%s: %g A Soll',field_name,local.excitation(g)));
        grid(heat_ax,'on');
        colormap(heat_ax,turbo(256));
        cb = colorbar(heat_ax);
        cb.Label.String = 'Faktor in A/A';
    end
end
