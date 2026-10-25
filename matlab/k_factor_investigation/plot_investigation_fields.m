function fig = plot_investigation_fields(a, export_format)
% The same field pipeline for every constant-shift method. Draw tabs on demand.
if nargin < 2
    export_format = '';
end
[fig,root_tab] = investigation_tab(a.set,'Feldrekonstruktion','Methoden');
method_tabs = uitabgroup(root_tab);
names = {'Messwerte','Zusammengeschoben','Zusammengeschoben + gemittelt', ...
    'Rekonstruktion','Messung + Rekonstruktion','Fehler','k-Verlauf'};
method_names = {'Kantenmittel','Längenmittel','Median','Trim 10 %','Bogen: Mittel', ...
    'Bogen: Median','Bogen: Trim','Segmentmittel','Dünn: Kanten','Dünn: Länge', ...
    'fixed iq','Diff. psi_d','Diff. psi_q','Diff. LS','Global F'};
ng = numel(a.excitation);
nm = numel(a.methods);
shown = true(ng,1);
heat_group = 1;
views = cell(nm,numel(names));
tabs = gobjects(nm,numel(names));
groups = a.groups;
x = a.x;
y = a.flux;
group_names = compose('%g A (Soll)',a.excitation);
components = {'d','q'};
scale = max(abs(x(:,1:2)),[],1);
triangles = cell(ng,1);
for g = 1:ng
    rows = groups==g;
    xy = x(rows,1:2)./scale;
    dt = delaunayTriangulation(xy);
    tri = dt.ConnectivityList;
    edges = dt.edges;
    limit = 9*mean(sum((xy(edges(:,1),:)-xy(edges(:,2),:)).^2,2));
    good = true(size(tri,1),1);
    for edge = [1 2;1 3;2 3]'
        good = good & sum((xy(tri(:,edge(1)),:)-xy(tri(:,edge(2)),:)).^2,2)<limit;
    end
    triangles{g} = tri(good,:);
end
error_limit = eps;
for m = 1:nm
    error_limit = max(error_limit,max(1000*abs(a.methods(m).error),[],'all','omitnan'));
end
error_colors = [linspace(.1,1,128)' linspace(.3,1,128)' ones(128,1); ...
    ones(128,1) linspace(1,.2,128)' linspace(1,.1,128)'];
ie_limits = [min(x(:,3)) max(x(:,3))];
rgb = interp1(linspace(ie_limits(1),ie_limits(2),256),turbo(256),x(:,3));
for m = 1:nm
    method_tab = uitab(method_tabs,'Title',method_names{m},'UserData',m, ...
        'Tag',sprintf('field_method_%02d',m));
    field_tabs = uitabgroup(method_tab);
    for s = 1:numel(names)
        tabs(m,s) = uitab(field_tabs,'Title',names{s},'UserData',s);
    end
    field_tabs.SelectionChangedFcn = @(~,event) build_view(m,event.NewValue.UserData);
end
method_tabs.SelectionChangedFcn = @(~,event) build_method(event.NewValue.UserData);
build_view(1,1);
if ~isempty(export_format)
    folder = fullfile(fileparts(mfilename('fullpath')),'results');
    for m = 1:nm
        for s = 1:numel(names)
            build_view(m,s);
            export_investigation_plot(views{m,s}.layout, ...
                fullfile(folder,sprintf('fields_set%d_method%02d_stage%02d',a.set,m,s)),export_format);
        end
    end
end

    function build_method(m)
        selected = tabs(m,1).Parent.SelectedTab.UserData;
        build_view(m,selected);
    end

    function build_view(m,s)
        if ~isempty(views{m,s})
            return;
        end
        if s==7
            layout = tiledlayout(tabs(m,s),1,2,'TileSpacing','compact');
            if m<=10
                values = a.ratios(:,:,m);
                label = 'abs(Delta id) / Delta ie';
            elseif m==11
                values = a.fixed_k;
                label = 'fixed iq: signed k, über gültige iq gemittelt';
            else
                values = a.methods(m).p(1)*ones(ng-1,numel(a.levels));
                label = 'Globaler konstanter d-Faktor';
            end
            excitation = (a.excitation(1:end-1)+a.excitation(2:end))/2;
            [L,E] = meshgrid(1000*a.levels,excitation);
            ax = nexttile(layout);
            surf(ax,L,E,values,'FaceColor','none','EdgeColor','interp');
            view(ax,35,25);
            zlabel(ax,'k in A/A');
            ax2 = nexttile(layout);
            imagesc(ax2,1000*a.levels,1:ng-1,values);
            yticks(ax2,1:ng-1);
            yticklabels(ax2,compose('%g -> %g A',a.excitation(1:end-1),a.excitation(2:end)));
            for p = [ax ax2]
                xlabel(p,'Konturlevel psi_d in mVs');
                grid(p,'on');
                colorbar(p);
                title(p,label);
            end
            ylabel(ax,'Mittlerer nominaler Erregerstrom des Paars in A');
            title(layout,[a.methods(m).name ': alle benachbarten Gruppen, 25–75 %']);
            views{m,s} = struct('layout',layout);
            return;
        end
        panel = uipanel(tabs(m,s),'BorderType','none','Position',[0 0 .85 1]);
        controls = uipanel(tabs(m,s),'Title','Erregerstromgruppen','Position',[.85 .02 .15 .96]);
        if s==6
            layout = tiledlayout(panel,2,3,'TileSpacing','compact');
        else
            layout = tiledlayout(panel,2,2,'TileSpacing','compact');
        end
        v = struct('layout',layout,'heat_axes',gobjects(2,1),'hist_axes',gobjects(0));
        if s~=3
            uicontrol(controls,'Style','pushbutton','String','Alle','Units','normalized', ...
                'Position',[.04 .92 .43 .05],'Callback',@(~,~) show_all(true));
            uicontrol(controls,'Style','pushbutton','String','Keine','Units','normalized', ...
                'Position',[.53 .92 .43 .05],'Callback',@(~,~) show_all(false));
            height = min(.045,.65/ng);
            for g = 1:ng
                uicontrol(controls,'Style','checkbox','String',group_names{g},'Value',shown(g), ...
                    'Tag',sprintf('reduction_check_%02d',g),'Units','normalized', ...
                    'Position',[.04 .87-g*height .94 height], ...
                    'Callback',@(source,~) show_group(g,logical(source.Value)));
            end
            uicontrol(controls,'Style','text','String','Heatmap-Gruppe','Units','normalized', ...
                'Position',[.04 .14 .92 .035]);
            uicontrol(controls,'Style','popupmenu','String',group_names,'Value',heat_group, ...
                'Tag','reduction_heatmap_group','Units','normalized','Position',[.04 .095 .92 .045], ...
                'Callback',@(source,~) select_heatmap(source.Value));
            note = {'Auswahl gilt in allen Methoden','und Feld-Tabs.','Nur Anzeige, keine Neuberechnung.'};
        else
            note = {'Gemeinsames Feld aus','allen Erregerstromgruppen.','Die Auswahl ändert nicht','das berechnete Mittelungsfeld.'};
        end
        uicontrol(controls,'Style','text','String',note,'Units','normalized','Position',[.02 .005 .96 .085]);
        for c = 1:2
            if s==6
                ax = nexttile(layout,(c-1)*3+1);
                v.heat_axes(c) = nexttile(layout,(c-1)*3+2);
                v.hist_axes(c) = nexttile(layout,(c-1)*3+3);
            else
                ax = nexttile(layout,(c-1)*2+1);
                v.heat_axes(c) = nexttile(layout,(c-1)*2+2);
            end
            hold(ax,'on');
            columns = [1 2];
            if c==1
                columns = [2 1]; % Rotate id/iq for the d-flux view.
            end
            if s==3
                map = a.methods(m).map;
                d = map.origin(1)+map.scale(1)*map.axes{1};
                q = map.origin(2)+map.scale(2)*map.axes{2};
                [D,Q] = ndgrid(d,q);
                coords = {D,Q};
                surf(ax,coords{columns(1)},coords{columns(2)},1000*map.values{c}, ...
                    'FaceColor','none','EdgeColor','interp');
                colormap(ax,turbo(256));
                clim(ax,1000*[min(y(:,c)) max(y(:,c))]);
                colorbar(ax);
            else
                coords = x(:,1:2);
                if s==2
                    coords = a.methods(m).coordinates;
                end
                values = y(:,c);
                if s==4
                    values = a.methods(m).prediction(:,c);
                elseif s==6
                    values = a.methods(m).error(:,c);
                end
                for g = 1:ng
                    rows = groups==g;
                    p = coords(rows,:);
                    z = 1000*values(rows);
                    tri = triangles{g};
                    tri = tri(all(isfinite(z(tri)),2),:);
                    color = x(rows,3);
                    if s==6
                        color = z;
                    end
                    h = trisurf(tri,p(:,columns(1)),p(:,columns(2)),z,color,'Parent',ax, ...
                        'FaceColor','none','EdgeColor','interp','LineWidth',.5);
                    h.Tag = sprintf('reduction_group_%02d',g);
                    if s==5
                        set(h,'FaceColor','interp','FaceVertexCData',rgb(rows,:), ...
                            'EdgeColor','none','FaceAlpha',.6);
                        predicted = 1000*a.methods(m).prediction(rows,c);
                        tri = tri(all(isfinite(predicted(tri)),2),:);
                        h2 = trisurf(tri,p(:,columns(1)),p(:,columns(2)),predicted,'Parent',ax, ...
                            'FaceColor','none','EdgeColor',[.12 .12 .12],'EdgeAlpha',.4,'LineWidth',.35);
                        h2.Tag = h.Tag;
                        h3 = scatter3(ax,p(:,columns(1)),p(:,columns(2)),predicted,18, ...
                            1000*a.methods(m).error(rows,c),'filled','MarkerEdgeColor','none','Tag',h.Tag);
                        if ~shown(g)
                            set([h2 h3],'Visible','off');
                        end
                    end
                    if ~shown(g)
                        h.Visible = 'off';
                    end
                end
                cb = colorbar(ax);
                if s==5
                    measured = plot3(ax,nan,nan,nan,'s','MarkerFaceColor',[.2 .5 .9], ...
                        'MarkerEdgeColor','none');
                    reconstructed = plot3(ax,nan,nan,nan,'-o','Color',[.12 .12 .12], ...
                        'MarkerFaceColor','w');
                    legend(ax,[measured reconstructed],{'Messung (Fläche)','Rekonstruktion (Fehlerfarbe)'}, ...
                        'Location','southoutside','Orientation','horizontal','FontSize',8,'AutoUpdate','off');
                end
                if ismember(s,[5 6])
                    colormap(ax,error_colors);
                    clim(ax,[-error_limit error_limit]);
                    cb.Label.String = 'Rekonstruktion - Messung in mVs';
                else
                    colormap(ax,turbo(256));
                    clim(ax,ie_limits);
                    cb.Label.String = 'Gemessener i_e in A';
                end
            end
            labels = coordinate_labels(m,s);
            xlabel(ax,labels{columns(1)});
            ylabel(ax,labels{columns(2)});
            zlabel(ax,sprintf('psi_%s in mVs',components{c}));
            if s==6
                zlabel(ax,'Rekonstruktion - Messung in mVs');
            end
            title(ax,sprintf('psi_%s: 3D',components{c}));
            view(ax,35,25);
            grid(ax,'on');
        end
        p = a.methods(m).p;
        title(layout,sprintf('%s | %s | konstant: [%s] A/A', ...
            a.methods(m).name,names{s},strtrim(sprintf('%.5g ',p))),'Interpreter','none');
        views{m,s} = v;
        draw_heatmaps(m,s);
    end

    function labels = coordinate_labels(m,s)
        labels = {'i_d in A','i_q in A'};
        if ismember(s,[2 3])
            labels{1} = 's_d = i_d + k_d i_e in A';
            if a.methods(m).model==6
                labels{2} = 's_q = i_q + k_q i_e in A';
            end
        end
    end

    function draw_heatmaps(m,s)
        v = views{m,s};
        g = heat_group;
        rows = find(groups==g);
        coords = x(rows,1:2);
        if s==2
            coords = a.methods(m).coordinates(rows,:);
        end
        labels = coordinate_labels(m,s);
        for c = 1:2
            ax = v.heat_axes(c);
            cla(ax);
            if s==3
                map = a.methods(m).map;
                d = map.origin(1)+map.scale(1)*map.axes{1};
                q = map.origin(2)+map.scale(2)*map.axes{2};
                [D,Q] = ndgrid(d,q);
                surf(ax,D,Q,zeros(size(D)),1000*map.values{c},'EdgeColor','none');
                title(ax,'Zusammengeschoben + gemittelt: gemeinsames Feld');
            else
                values = y(rows,c);
                if s==4
                    values = a.methods(m).prediction(rows,c);
                elseif ismember(s,[5 6])
                    values = a.methods(m).error(rows,c);
                end
                values = 1000*values;
                if ~shown(g)
                    values(:) = nan;
                end
                tri = triangles{g};
                tri = tri(all(isfinite(values(tri)),2),:);
                trisurf(tri,coords(:,1),coords(:,2),zeros(size(values)),values,'Parent',ax, ...
                    'FaceColor','interp','EdgeColor','none');
                title(ax,sprintf('Heatmap: %g A (Soll)',a.excitation(g)));
                if s==6
                    h = v.hist_axes(c);
                    cla(h);
                    samples = values(isfinite(values));
                    histogram(h,samples,'Normalization','pdf','FaceColor',[.2 .45 .75]);
                    hold(h,'on');
                    if ~isempty(samples)
                        xline(h,mean(samples),'k-');
                        xline(h,mean(samples)-std(samples),'k--');
                        xline(h,mean(samples)+std(samples),'k--');
                    end
                    title(h,{sprintf('Mittel %.3g, Std %.3g mVs',mean(samples),std(samples)), ...
                        sprintf('RMSE %.3g mVs, n=%d',sqrt(mean(samples.^2)),numel(samples))});
                    xlabel(h,'Fehler in mVs');
                    ylabel(h,'Wahrscheinlichkeitsdichte');
                    grid(h,'on');
                end
            end
            if ismember(s,[5 6])
                colormap(ax,error_colors);
                clim(ax,[-error_limit error_limit]);
            else
                colormap(ax,turbo(256));
                clim(ax,1000*[min(y(:,c)) max(y(:,c))]);
            end
            cb = colorbar(ax);
            if ismember(s,[5 6])
                cb.Label.String = 'Rekonstruktion - Messung in mVs';
            else
                cb.Label.String = ['psi_' components{c} ' in mVs'];
            end
            view(ax,2);
            axis(ax,'tight');
            xlabel(ax,labels{1});
            ylabel(ax,labels{2});
        end
    end

    function show_group(g,value)
        shown(g) = value;
        state = 'off';
        if value
            state = 'on';
            heat_group = g;
        end
        set(findall(root_tab,'Tag',sprintf('reduction_group_%02d',g)),'Visible',state);
        set(findall(root_tab,'Tag',sprintf('reduction_check_%02d',g)),'Value',double(value));
        refresh_heatmaps;
    end

    function show_all(value)
        shown(:) = value;
        state = 'off';
        if value
            state = 'on';
        end
        for g = 1:ng
            set(findall(root_tab,'Tag',sprintf('reduction_group_%02d',g)),'Visible',state);
            set(findall(root_tab,'Tag',sprintf('reduction_check_%02d',g)),'Value',double(value));
        end
        refresh_heatmaps;
    end

    function select_heatmap(g)
        heat_group = g;
        show_group(g,true);
    end

    function refresh_heatmaps
        set(findall(root_tab,'Tag','reduction_heatmap_group'),'Value',heat_group);
        for m = 1:nm
            for s = [1 2 4 5 6]
                if ~isempty(views{m,s})
                    draw_heatmaps(m,s);
                end
            end
        end
    end
end
