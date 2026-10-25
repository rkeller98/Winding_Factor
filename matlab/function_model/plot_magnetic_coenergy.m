function [fig,data] = plot_magnetic_coenergy(dataset)
% Interactive saved coenergy/flux/inductance fields; no fits or paper exports.
% addpath('C:\Git\Winding_Factor\matlab\function_model');
% plot_magnetic_coenergy             % Set2, all excitation groups
% plot_magnetic_coenergy('Set1')     % separate recording
% plot_magnetic_coenergy('pooled')   % existing joint fit, groups kept separate
% W', psi_e and L_ee use W'(0,0,ie)=0; their absolute values are unidentified.
if nargin < 1
    dataset = 'Set2';
end
data = prepare_magnetic_coenergy(dataset);
names = {'W'' rel','Psi_d','Psi_q','Psi_e rel', ...
    'L_dd','L_dq','L_de','L_qd','L_qq','L_qe','L_ed','L_eq','L_ee rel'};
symbols = {'W''_{rel}','\Psi_d','\Psi_q','\Psi_{e,rel}', ...
    'L_{dd}','L_{dq}','L_{de}','L_{qd}','L_{qq}','L_{qe}', ...
    'L_{ed}','L_{eq}','L_{ee,rel}'};
units = [{'J'},repmat({'mVs'},1,3),repmat({'mH'},1,9)];
factor = [1,repmat(1000,1,12)];
ng = numel(data.labels);
shown = true(ng,1);
method = 3;
mode = 1;
heat_group = 1;
current = 1;
group_colors = lines(ng);
fig = figure('Name',['Magnetische Koenergie | ' data.name], ...
    'NumberTitle','off','Position',[40 40 1500 900],'Color','w');
theme(fig,'light');
tabs = uitabgroup(fig,'Units','normalized','Position',[.01 .17 .78 .81]);
tab = gobjects(13,1);
for k = 1:13
    tab(k) = uitab(tabs,'Title',names{k},'UserData',k);
end
axes3 = gobjects(13,1);
axes2 = gobjects(13,1);
keys = repmat({''},13,1);
panel = uipanel(fig,'Title','Auswahl (alle Tabs)','Units','normalized', ...
    'Position',[.80 .02 .19 .96]);
uicontrol(panel,'Style','text','String','Gespeichertes Forschungsmodell', ...
    'Units','normalized','Position',[.03 .94 .94 .025]);
uicontrol(panel,'Style','popupmenu','String',data.methods,'Value',method, ...
    'Units','normalized','Position',[.03 .90 .94 .04], ...
    'Tag','coenergy_method','Callback',@change_method);
uicontrol(panel,'Style','text','String','Psi_d / Psi_q: Datenansicht', ...
    'Units','normalized','Position',[.03 .87 .94 .025]);
uicontrol(panel,'Style','popupmenu', ...
    'String',{'Modell','Messwerte (GUI-ausgewertet)','Modell - Messwerte', ...
    'Messung: q-Spiegelabweichung'}, ...
    'Units','normalized','Position',[.03 .83 .94 .04], ...
    'Tag','coenergy_mode','Callback',@change_mode);
uicontrol(panel,'Style','pushbutton','String','Alle', ...
    'Units','normalized','Position',[.03 .78 .45 .04], ...
    'Tag','coenergy_all','Callback',@(~,~) show_all(true));
uicontrol(panel,'Style','pushbutton','String','Keine', ...
    'Units','normalized','Position',[.52 .78 .45 .04], ...
    'Tag','coenergy_none','Callback',@(~,~) show_all(false));
groups_panel = uipanel(panel,'BorderType','none','Units','normalized', ...
    'Position',[.01 .24 .98 .53]);
checks = gobjects(ng,1);
for g = 1:ng
    checks(g) = uicontrol(groups_panel,'Style','checkbox', ...
        'String',data.labels{g},'Value',1,'FontSize',8, ...
        'Units','normalized','Position',[.02 1-g/ng .97 1/ng], ...
        'ForegroundColor',.7*group_colors(g,:), ...
        'Tag',sprintf('coenergy_check_%d',g), ...
        'Callback',@(source,~) show_group(g,logical(source.Value)));
end
uicontrol(panel,'Style','text','String','Heatmap: eine sichtbare Gruppe', ...
    'Units','normalized','Position',[.03 .20 .94 .03]);
heat_popup = uicontrol(panel,'Style','popupmenu','String',data.labels, ...
    'Units','normalized','Position',[.03 .16 .94 .04], ...
    'Tag','coenergy_heat_group','Callback',@change_heat);
uicontrol(panel,'Style','pushbutton','String','Aktuellen Tab als PNG speichern', ...
    'Units','normalized','Position',[.03 .10 .94 .045], ...
    'Tag','coenergy_export','Callback',@save_png);
uicontrol(panel,'Style','text','String', ...
    {'3D: drehen / zoomen mit Maus.','Farbskala je Groesse fest fuer alle', ...
     'Gruppen und vier Modelle.'}, ...
    'Units','normalized','Position',[.03 .01 .94 .08]);
info = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[.02 .01 .76 .15],'HorizontalAlignment','left', ...
    'BackgroundColor','w','FontSize',10,'Tag','coenergy_info');
tabs.SelectionChangedFcn = @(~,event) select_tab(event.NewValue.UserData);
rotate3d(fig,'on');
draw_tab();
fprintf('%s: %d points, %d separate excitation groups. No refit.\n', ...
    data.name,size(data.x,1),ng);
fprintf('Source: %s\n',data.source);
fprintf('W'', Psi_e, L_ee: relative; %d points hidden (reference outside measured hull).\n', ...
    sum(~data.reference_supported));

    function select_tab(k)
        current = k;
        draw_tab();
    end

    function change_method(source,~)
        method = source.Value;
        draw_tab();
    end

    function change_mode(source,~)
        mode = source.Value;
        draw_tab();
    end

    function change_heat(source,~)
        heat_group = source.Value;
        % A hidden group stays hidden in both views.
        draw_tab();
    end

    function show_all(state)
        shown(:) = state;
        set(checks,'Value',double(state));
        update_visibility();
    end

    function show_group(g,state)
        shown(g) = state;
        update_visibility();
    end

    function update_visibility()
        if ~shown(heat_group) && any(shown)
            heat_group = find(shown,1);
            heat_popup.Value = heat_group;
        end
        for group = 1:ng
            objects = findall(tabs,'Tag',sprintf('coenergy_group_%d',group));
            set(objects,'Visible',onoff(shown(group)));
        end
        draw_tab();
    end

    function [z,limits,description] = values(k)
        z = data.values{method}(:,k)*factor(k);
        all_values = cellfun(@(v) v(:,k)*factor(k),data.values,'UniformOutput',false);
        description = [data.methods{method} ': rekonstruierter Modellwert'];
        if k == 2 || k == 3
            measured = data.measured(:,k-1)*factor(k);
            switch mode
                case 1
                    all_values{end+1} = measured;
                    description = [description '; schwarze Punkte: Messwerte'];
                case 2
                    z = measured;
                    all_values = {measured};
                    description = 'GUI-ausgewertete Statorfluss-Messwerte, unveraendert';
                case 3
                    z = z-measured;
                    all_values = cellfun(@(v) v-measured,all_values,'UniformOutput',false);
                    description = 'Modell minus Messwert am selben gemessenen (id,iq,ie)';
                case 4
                    z = 1000*data.mirror_residual(:,k-1);
                    all_values = {z};
                    description = ['Messsymmetrie: Psi_d(id,iq)-Psi_d(id,-iq); ' ...
                        'Psi_q(id,iq)+Psi_q(id,-iq). Keine Symmetrie erzwungen.'];
            end
        end
        joined = vertcat(all_values{:});
        finite = joined(isfinite(joined));
        if isempty(finite)
            limits = [-1 1];
        elseif mode >= 3 && (k == 2 || k == 3)
            amplitude = max(max(abs(finite)),eps);
            limits = [-amplitude amplitude];
        else
            limits = [min(finite),max(finite)];
            if limits(1) == limits(2)
                limits = limits+[-1 1]*max(abs(limits(1))*.01,1e-6);
            end
        end
    end

    function draw_tab()
        k = current;
        [z,limits,description] = values(k);
        symbol = symbols{k};
        if (k == 2 || k == 3) && mode >= 3
            symbol = ['\Delta ' symbol];
        end
        key = sprintf('%d/%d/%d',method,mode,heat_group);
        if ~strcmp(keys{k},key)
            delete(tab(k).Children);
            axes3(k) = axes(tab(k),'Position',[.07 .14 .37 .74]);
            axes2(k) = axes(tab(k),'Position',[.57 .14 .34 .74]);
            hold(axes3(k),'on');
            for group = 1:ng
                rows = find(data.groups == group);
                tri = data.triangles{group};
                tri = tri(all(isfinite(z(tri)),2),:);
                tag = sprintf('coenergy_group_%d',group);
                if ~isempty(tri)
                    patch(axes3(k),'Faces',tri,'Vertices',[data.x(:,1:2),z], ...
                        'FaceColor','none','EdgeColor',group_colors(group,:), ...
                        'LineWidth',.6,'Tag',tag,'Visible',onoff(shown(group)));
                end
                good = rows(isfinite(z(rows)));
                scatter3(axes3(k),data.x(good,1),data.x(good,2),z(good), ...
                    7,z(good),'filled','Tag',tag,'Visible',onoff(shown(group)));
                if (k == 2 || k == 3) && mode == 1
                    scatter3(axes3(k),data.x(rows,1),data.x(rows,2), ...
                        1000*data.measured(rows,k-1),4,[.15 .15 .15], ...
                        'Tag',tag,'Visible',onoff(shown(group)));
                end
            end
            group = heat_group;
            rows = find(data.groups == group);
            tri = data.triangles{group};
            tri = tri(all(isfinite(z(tri)),2),:);
            tag = sprintf('coenergy_group_%d',group);
            if ~isempty(tri)
                patch(axes2(k),'Faces',tri,'Vertices', ...
                    [data.x(:,1:2),zeros(size(z))],'FaceVertexCData',z, ...
                    'FaceColor','interp','EdgeColor','none', ...
                    'Tag',tag,'Visible',onoff(shown(group)));
            end
            hold(axes2(k),'on');
            good = rows(isfinite(z(rows)));
            scatter(axes2(k),data.x(good,1),data.x(good,2),8,z(good), ...
                'filled','Tag',tag,'Visible',onoff(shown(group)));
            for ax = [axes3(k),axes2(k)]
                xlabel(ax,'i_d in A');
                ylabel(ax,'i_q in A');
                grid(ax,'on');
                colormap(ax,turbo(256));
                clim(ax,limits);
                cb = colorbar(ax);
                cb.Label.String = [symbol ' in ' units{k}];
                xlim(ax,[min(data.x(:,1)),max(data.x(:,1))]);
                ylim(ax,[min(data.x(:,2)),max(data.x(:,2))]);
            end
            zlabel(axes3(k),[symbol ' in ' units{k}]);
            zlim(axes3(k),limits);
            view(axes3(k),35,25);
            view(axes2(k),2);
            title(axes3(k),{[data.name ' | ' symbol], ...
                'Meshfarbe: Gruppe; Punktfarbe: Wert'},'FontSize',11);
            title(axes2(k),{symbol,data.labels{heat_group}}, ...
                'Interpreter','tex','FontSize',10);
            keys{k} = key;
        end
        update_info(description);
    end

    function update_info(description)
        relative = '';
        if ismember(current,data.relative)
            relative = sprintf(['RELATIV: W''(0,0,ie)=0; C(ie), C''(ie), C''''(ie) unbekannt. ' ...
                '%d Punkte mit Referenz ausserhalb der Messhuelle ausgeblendet.\n'], ...
                sum(~data.reference_supported));
        end
        if ~any(shown)
            state = 'Keine Erregerstromgruppe sichtbar.';
        elseif ~shown(heat_group)
            state = 'Heatmap-Gruppe ausgeblendet.';
        else
            state = ['Heatmap: ' data.labels{heat_group}];
        end
        if (current == 2 || current == 3) && mode == 4
            delta_ie = data.mirror_ie_difference(data.groups == heat_group);
            delta_ie = delta_ie(isfinite(delta_ie));
            state = sprintf(['%s | Spiegelpaare derselben Sollgruppe, ie kann abweichen: ' ...
                'RMS Delta ie = %.3g A (%d Punkte).'],state, ...
                sqrt(mean(delta_ie.^2)),numel(delta_ie));
        end
        info.String = sprintf(['%s\n%s\n%s' ...
            'Alle Flaechen an tatsaechlichen Messstroemen, kein konstanter ie-Schnitt. ' ...
            'Heatmap: lineare Darstellung innerhalb getrennter Gruppen-Dreiecke.\n' ...
            'Alle L sind Modellableitungen, keine unabhaengigen L-Messungen. ' ...
            'D L symmetrisch: Ldq=Lqd, Led=1.5 Lde, Leq=1.5 Lqe.\n' ...
            'Keine q-Spiegelung erzwungen. Reziprozitaet des Potentials belegt keine Messsymmetrie.'], ...
            description,state,relative);
    end

    function save_png(~,~)
        [file,folder] = uiputfile('*.png','Aktuelle Ansicht speichern', ...
            [strrep(data.name,' + ','_') '_' data.methods{method} '_' ...
            regexprep(names{current},'[^A-Za-z0-9]','') '.png']);
        if isequal(file,0)
            return
        end
        % Export the current tab plus the visible physical interpretation.
        exportapp(fig,fullfile(folder,file));
    end
end

function state = onoff(value)
if value
    state = 'on';
else
    state = 'off';
end
end
