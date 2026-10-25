function overview = plot_investigation_comparison(study, export_format)
% Compare each method within its own criterion; no new fit or flux prediction.
if nargin < 2
    export_format = '';
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab','common'));
folder = fullfile(root,'matlab','k_factor_investigation','results');
nsets = numel(study);
set_ids = cellfun(@(a) a.set,study);
global_factors = nan(nsets,5);
residuals = nan(nsets,5);
pair_factors = nan(nsets,5);
factor_names = {'k aus psi_d','k aus psi_q','lokales LS','global Kontur','global Differential'};
for j = 1:nsets
    item = study{j};
    c = item.contour;
    l = item.local;
    table = l.component_comparison;
    global_factors(j,:) = table.median';
    residuals(j,:) = table.relative_derivative_residual';
    pair_factors(j,:) = [c.signed_k(1:4)' median(c.fixed_iq.signed_k,'omitnan')];

    [~,tab] = investigation_tab(item.set,'Methodenvergleich','Konturverfahren');
    layout = tiledlayout(tab,1,2,'TileSpacing','compact');
    ax = nexttile(layout);
    k = [c.signed_k;median(c.fixed_iq.signed_k,'omitnan')];
    names = [c.methods {'fixed iq (Median über iq)'}];
    plot(ax,k,1:numel(k),'o','LineWidth',1.2);
    hold(ax,'on');
    fixed = c.fixed_iq.signed_k;
    fixed = fixed(isfinite(fixed));
    if ~isempty(fixed)
        plot(ax,[min(fixed) max(fixed)],[numel(k) numel(k)],'k-','LineWidth',1.5);
    end
    xline(ax,l.global_contour_k,'k--','globaler Konturfaktor');
    yticks(ax,1:numel(k));
    yticklabels(ax,names);
    ylim(ax,[.5 numel(k)+.5]);
    xlabel(ax,'Vorzeichenbehaftetes k in A/A');
    grid(ax,'on');
    title(ax,{'Ein Konturlevel / ein Gruppenpaar', ...
        'fixed iq: Strich zeigt Min/Max über gültige iq'});
    ax = nexttile(layout);
    imagesc(ax,1000*c.level_comparison.levels,1:numel(c.methods), ...
        c.level_comparison.signed_k');
    yticks(ax,1:numel(c.methods));
    yticklabels(ax,c.methods);
    xticks(ax,1000*c.level_comparison.levels);
    xlabel(ax,'Konturlevel psi_d in mVs');
    cb = colorbar(ax);
    cb.Label.String = 'Vorzeichenbehaftetes k in A/A';
    colormap(ax,turbo(256));
    title(ax,'Alle Mittelungsverfahren auf drei Flussniveaus');
    title(layout,sprintf('Set%d: nominal %g/%g A, gewähltes Niveau %.4g mVs; gemessenes Delta ie', ...
        item.set,c.nominal_excitation,1000*c.level));
    export_investigation_plot(layout,fullfile(folder,sprintf('method_contours_set%d',item.set)),export_format);

    [~,tab] = investigation_tab(item.set,'Methodenvergleich','Differentialfaktoren');
    layout = tiledlayout(tab,1,2,'TileSpacing','compact');
    ax = nexttile(layout);
    errorbar(ax,1:3,table.median(1:3),table.median(1:3)-table.p05(1:3), ...
        table.p95(1:3)-table.median(1:3),'o','LineWidth',1.4);
    hold(ax,'on');
    yline(ax,l.global_contour_k,'k--','Kontur global');
    yline(ax,l.global_derivative_k,'k:','Differential global');
    xticks(ax,1:3);
    xticklabels(ax,factor_names(1:3));
    ylabel(ax,'k in A/A');
    grid(ax,'on');
    title(ax,{'5/50/95%-Quantile auf identischen Punkten', ...
        'Kennfeldstreuung, keine Unsicherheitsintervalle'});
    ax = nexttile(layout);
    bar(ax,residuals(j,:));
    xticks(ax,1:5);
    xticklabels(ax,factor_names);
    xtickangle(ax,25);
    ylabel(ax,'||b - k a||_F / ||b||_F');
    grid(ax,'on');
    title(ax,{'Gleiches Differentialkriterium / gleiche Punkte', ...
        'Innerhalb der Ableitungsdaten, kein Flussvorhersagefehler'});
    title(layout,sprintf('Set%d: %d gemeinsame gültige Punkte; a=[Ldd,Lqd], b=[Lde,Lqe]', ...
        item.set,nnz(l.common_scalar)));
    export_investigation_plot(layout,fullfile(folder,sprintf('method_differential_set%d',item.set)),export_format);
    if isfield(item,'fields')
        f = item.fields;
        [~,tab] = investigation_tab(item.set,'Methodenvergleich','Feldrekonstruktion');
        layout = tiledlayout(tab,1,2,'TileSpacing','compact');
        ax = nexttile(layout);
        barh(ax,f.rmse);
        yticks(ax,1:numel(f.methods));
        yticklabels(ax,{f.methods.name});
        xlabel(ax,'RMSE auf gemeinsamer Abdeckung in mVs');
        legend(ax,'psi_d','psi_q','Location','best');
        grid(ax,'on');
        title(ax,sprintf('%d identische gültige Messpunkte',nnz(f.common)));
        ax = nexttile(layout);
        coverage = arrayfun(@(m) nnz(m.valid),f.methods)/numel(f.groups);
        barh(ax,100*coverage);
        yticks(ax,1:numel(f.methods));
        yticklabels(ax,{f.methods.name});
        xlabel(ax,'Gültig rekonstruierte Messpunkte in %');
        xline(ax,100*mean(f.common),'k--','gemeinsame Abdeckung');
        grid(ax,'on');
        title(ax,'Abdeckung getrennt vom Fehler vergleichen');
        title(layout,'Gleiche 81x81-LUT und Interpolation; Fehler innerhalb der verwendeten Messdaten');
        export_investigation_plot(layout,fullfile(folder,sprintf('field_errors_set%d',item.set)),export_format);
        for c = 1:2
            component = {'d','q'};
            [~,tab] = investigation_tab(item.set,'Methodenvergleich',['Gruppenfehler psi_' component{c}]);
            layout = tiledlayout(tab,1,2,'TileSpacing','compact');
            [M,E] = meshgrid(1:numel(f.methods),f.excitation);
            values = f.group_rmse(:,:,c);
            ax = nexttile(layout);
            surf(ax,M,E,values,'FaceColor','none','EdgeColor','interp');
            view(ax,35,25);
            zlabel(ax,'RMSE in mVs');
            ax2 = nexttile(layout);
            imagesc(ax2,1:numel(f.methods),1:numel(f.excitation),values);
            yticks(ax2,1:numel(f.excitation));
            yticklabels(ax2,compose('%g A',f.excitation));
            for p = [ax ax2]
                xticks(p,1:numel(f.methods));
                xticklabels(p,{f.methods.name});
                xtickangle(p,45);
                ylabel(p,'Nominaler Erregerstrom');
                colorbar(p);
                colormap(p,turbo(256));
                clim(p,[0 max(max(values,[],'all','omitnan'),eps)]);
            end
            title(layout,sprintf('psi_%s: alle Gruppen auf gemeinsamer Auswertemaske',component{c}));
            export_investigation_plot(layout,fullfile(folder,sprintf('field_group_errors_set%d_%s',item.set,component{c})),export_format);
        end
    end
end

[overview,tab] = investigation_tab('Methodenvergleich','Alle Sets','Faktoren und Differentialkriterium');
layout = tiledlayout(tab,1,2,'TileSpacing','compact');
ax = nexttile(layout);
plot(ax,1:5,global_factors','o-','LineWidth',1.2);
xticks(ax,1:5);
xticklabels(ax,factor_names);
xtickangle(ax,20);
ylabel(ax,'k in A/A');
grid(ax,'on');
legend(ax,compose('Set%d',set_ids),'Location','best');
title(ax,'Lokale Mediane / globale Faktoren je Set');
ax = nexttile(layout);
bar(ax,residuals');
xticks(ax,1:5);
xticklabels(ax,factor_names);
xtickangle(ax,20);
ylabel(ax,'Relativer Differentialrest');
grid(ax,'on');
legend(ax,compose('Set%d',set_ids),'Location','best');
title(ax,'Gleiche Punkte je Set, unterschiedliche Stützbereiche zwischen Sets');
title(layout,'Methoden nach demselben Differentialkriterium; keine neue Optimierung');
export_investigation_plot(layout,fullfile(folder,'method_overview_differential'),export_format);

[~,tab] = investigation_tab('Methodenvergleich','Alle Sets','Konturen und Widerstand');
layout = tiledlayout(tab,1,2,'TileSpacing','compact');
ax = nexttile(layout);
plot(ax,1:5,pair_factors','o-','LineWidth',1.2);
xticks(ax,1:5);
xticklabels(ax,[study{1}.contour.methods(1:4) {'fixed iq (Median)'}]);
xtickangle(ax,20);
ylabel(ax,'Vorzeichenbehaftetes k in A/A');
grid(ax,'on');
legend(ax,compose('Set%d',set_ids),'Location','best');
title(ax,{'Repräsentatives Niveau / Paar je Set', ...
    'Niveaus und Gruppen sind zwischen Sets unterschiedlich'});
ax = nexttile(layout);
hold(ax,'on');
for j = 1:nsets
    r = study{j}.resistance;
    if r.skipped
        continue;
    end
    baseline = r.summary(r.offset_fractions==0,2);
    plot(ax,100*r.offset_fractions,100*(r.summary(:,2)/baseline-1),'o-', ...
        'DisplayName',sprintf('Set%d',study{j}.set));
end
xlabel(ax,'Rs-Offset in % des ursprünglichen medianen Rs');
ylabel(ax,'Änderung des globalen Konturfaktors in %');
grid(ax,'on');
legend(ax,'Location','best');
title(ax,'Widerstandssensitivität, keine Rs-Identifikation');
title(layout,'Konturverfahren separat vom Differentialkriterium vergleichen');
export_investigation_plot(layout,fullfile(folder,'method_overview_contours'),export_format);
if all(cellfun(@(s) isfield(s,'fields'),study))
    [~,tab] = investigation_tab('Methodenvergleich','Alle Sets','Feldrekonstruktion');
    layout = tiledlayout(tab,1,2,'TileSpacing','compact');
    for c = 1:2
        ax = nexttile(layout);
        values = nan(numel(study{1}.fields.methods),nsets);
        for j = 1:nsets
            values(:,j) = study{j}.fields.rmse(:,c);
        end
        barh(ax,values);
        yticks(ax,1:numel(study{1}.fields.methods));
        yticklabels(ax,{study{1}.fields.methods.name});
        xlabel(ax,'RMSE in mVs');
        legend(ax,compose('Set%d',set_ids),'Location','best');
        grid(ax,'on');
        component = {'psi_d','psi_q'};
        title(ax,component{c},'Interpreter','none');
    end
    title(layout,'Feldrekonstruktion: gleiche Punkte innerhalb jedes Sets; unterschiedliche Messbereiche zwischen Sets');
    export_investigation_plot(layout,fullfile(folder,'field_error_overview'),export_format);
end
end
