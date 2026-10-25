function figures = plot_local_winding_factors
% Display the saved joint local differential factor, without recalculating it.
folder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(folder,'common'));
root = fileparts(folder);
loaded = load(fullfile(root,'results','reduction_data.mat'),'data');
saved = load(fullfile(root,'results','winding_validation.mat'),'study');
figures = gobjects(numel(loaded.data),1);
for s = 1:numel(loaded.data)
    [x,~,groups,excitation] = reduction_inputs(loaded.data{s});
    figures(s) = plot_set(x,groups,excitation,saved.study{s}.slopes{2},s);
end
end

function fig = plot_set(x,groups,excitation,slopes,set_index)
% Mesh and group controls taken from k_factor_investigation/plot_local_factor.
values = slopes.klocal;
valid = slopes.valid & isfinite(values) & all(isfinite(x),2);
scale = max(abs(x(:,1:2)),[],1);
scale(scale==0) = 1;
ng = numel(excitation);
fig = figure('Visible','on','Color','w','Tag','local_winding_factor', ...
    'Name',sprintf('Set%d: local differential k-factor',set_index), ...
    'Position',[40 40 1400 850]);
theme(fig,'light');
ax = axes(fig,'Position',[.07 .12 .64 .78]);
hold(ax,'on');
meshes = gobjects(ng,1);
points = gobjects(ng,1);
checks = gobjects(ng,1);
panel = uipanel(fig,'Title','Erregerstromgruppen','Units','normalized', ...
    'Position',[.82 .04 .17 .92]);
uicontrol(panel,'Style','pushbutton','String','Alle','Tag','show_all_groups', ...
    'Units','normalized','Position',[.04 .92 .43 .05],'Callback',@(~,~) show_all(true));
uicontrol(panel,'Style','pushbutton','String','Keine','Tag','hide_all_groups', ...
    'Units','normalized','Position',[.53 .92 .43 .05],'Callback',@(~,~) show_all(false));
height = min(.05,.78/ng);
for g = 1:ng
    rows = find(groups==g);
    xy = x(rows,1:2)./scale;
    % Keep the original group mesh; invalid vertices leave holes.
    dt = delaunayTriangulation(xy);
    tri = dt.ConnectivityList;
    edges = dt.edges;
    limit = 9*mean(sum((xy(edges(:,1),:)-xy(edges(:,2),:)).^2,2));
    keep = all(valid(rows(tri)),2);
    for edge = [1 2;1 3;2 3]'
        keep = keep & sum((xy(tri(:,edge(1)),:)-xy(tri(:,edge(2)),:)).^2,2)<limit;
    end
    z = values(rows);
    z(~valid(rows)) = nan;
    meshes(g) = patch(ax,'Faces',tri(keep,:),'Vertices',[x(rows,1:2) z], ...
        'FaceVertexCData',x(rows,3),'FaceColor','none','EdgeColor','interp', ...
        'LineWidth',.5,'Tag',sprintf('local_mesh_%02d',g));
    rows = rows(valid(rows));
    points(g) = scatter3(ax,x(rows,1),x(rows,2),values(rows),9,x(rows,3), ...
        'filled','Tag',sprintf('local_points_%02d',g));
    checks(g) = uicontrol(panel,'Style','checkbox','String',sprintf('%g A (Soll)',excitation(g)), ...
        'Value',1,'Units','normalized','Position',[.04 .87-g*height .94 height], ...
        'Tag',sprintf('local_group_%02d',g), ...
        'Callback',@(source,~) show_group(g,logical(source.Value)));
end
uicontrol(panel,'Style','text','String',{'Auswahl ändert nur','die Sichtbarkeit.'}, ...
    'Units','normalized','Position',[.04 .01 .92 .07]);
xlabel(ax,'i_d in A');
ylabel(ax,'i_q in A');
zlabel(ax,'k_{loc} in A/A');
title(ax,{sprintf('Set%d: gemeinsamer lokaler Differentialfaktor (80 Nachbarn)',set_index), ...
    'k_{loc} = (L_{dd}L_{de} + L_{qd}L_{qe}) / (L_{dd}^2 + L_{qd}^2)'});
grid(ax,'on');
view(ax,35,25);
colormap(ax,turbo(256));
clim(ax,[min(x(:,3)) max(x(:,3))]);
cb = colorbar(ax);
cb.Position = [.75 .15 .016 .68];
cb.Label.String = 'Gemessener i_e in A';
rotate3d(fig,'on');

    function show_group(g,visible)
        state = 'off';
        if visible
            state = 'on';
        end
        set([meshes(g);points(g)],'Visible',state);
    end

    function show_all(visible)
        set(checks,'Value',double(visible));
        state = 'off';
        if visible
            state = 'on';
        end
        set([meshes;points],'Visible',state);
    end
end
