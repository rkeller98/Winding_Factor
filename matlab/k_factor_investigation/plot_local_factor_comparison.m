function fig = plot_local_factor_comparison(local, group, export_format)
% All excitation groups are selectable; spatial maps have both projections.
if nargin < 2 || isempty(group)
    group = 1;
end
if nargin < 3
    export_format = '';
end
[fig,plot_tab] = investigation_tab(local.set,'Lokale Faktoren','Differenzen und Gruppen');
selection = uitabgroup(plot_tab,'Position',[0 .94 1 .06]);
for g = 1:numel(local.excitation)
    tab = uitab(selection,'Title',sprintf('%g A',local.excitation(g)),'UserData',g);
    if g==group
        selection.SelectedTab = tab;
    end
end
selection.SelectionChangedFcn = @(~,event) plot_local_factor_comparison(local,event.NewValue.UserData,'');
plot_tabs = uitabgroup(plot_tab,'Position',[0 0 1 .94]);
x = local.x;
rows = local.groups==group;
masks = {rows & local.common_scalar,rows & isfinite(local.L(:,2,1)),rows & local.valid_local};
values = {local.kd_from_psi_q-local.kd_from_psi_d,1000*local.L(:,2,1),1000*local.L(:,1,1)};
names = {'k(psi_q) - k(psi_d) in A/A','Lqd in mH','Ldd in mH'};
layouts = gobjects(4,1);
for metric = 1:3
    tab = uitab(plot_tabs,'Title',names{metric});
    layouts(metric) = tiledlayout(tab,1,2,'TileSpacing','compact');
    for projection = 1:2
        ax = nexttile(layouts(metric));
        keep = masks{metric};
        v = values{metric};
        scatter3(ax,x(keep,1),x(keep,2),v(keep),18,v(keep),'filled');
        hold(ax,'on');
        if metric==2
            masked = rows & local.small_Lqd;
            scatter3(ax,x(masked,1),x(masked,2),v(masked),40,'k','o');
        end
        if projection==1
            view(ax,35,25);
        else
            view(ax,2);
        end
        xlabel(ax,'i_d in A');
        ylabel(ax,'i_q in A');
        zlabel(ax,names{metric},'Interpreter','none');
        title(ax,names{metric},'Interpreter','none');
        colorbar(ax);
        grid(ax,'on');
        colormap(ax,turbo(256));
    end
    title(layouts(metric),sprintf('Set%d, %g A Soll; schwarze Kreise: kleine Lqd, keine Wertbegrenzung', ...
        local.set,local.excitation(group)));
end
tab = uitab(plot_tabs,'Title','Alle Gruppen');
layouts(4) = tiledlayout(tab,1,1,'TileSpacing','compact');
ax = nexttile(layouts(4));
g = local.component_group_summary;
plot(ax,g(:,1),g(:,4:6),'.-','LineWidth',1.2);
hold(ax,'on');
yline(ax,local.global_contour_k,'--k');
yline(ax,local.global_derivative_k,':k');
xlabel(ax,'Nominaler i_e in A');
ylabel(ax,'Faktor in A/A');
grid(ax,'on');
legend(ax,{'k_d^{(\psi_d)}','k_d^{(\psi_q)}','k_{LS}','Global Kontur','Global Differential'},'Location','best');
title(ax,'Gruppenmediane der drei lokalen Faktoren auf identischen Punkten');
if ~isempty(export_format)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root,'matlab','common'));
    for k = 1:numel(layouts)
        export_investigation_plot(layouts(k),fullfile(root,'matlab','k_factor_investigation','results', ...
            sprintf('local_comparison_set%d_group%02d_%d',local.set,group,k)),export_format);
    end
end
end
