function fig = plot_local_factor_comparison(local, group, export_format)
% Targeted maps and group summaries; associations do not identify a cause.
if nargin < 2 || isempty(group)
    group = ceil(numel(local.excitation)/2);
end
if nargin < 3
    export_format = '';
end
rows = local.groups == group;
common = rows & local.common_scalar;
x = local.x;
[fig,plot_tab] = investigation_tab(local.set,'Lokale Faktoren','Differenzen und Gruppen');
t = tiledlayout(plot_tab, 2, 2, 'TileSpacing', 'compact');
ax = nexttile(t);
scatter(ax, x(common, 1), x(common, 2), 20, ...
    local.kd_from_psi_q(common)-local.kd_from_psi_d(common), 'filled');
title(ax, 'k_d^{(\psi_q)} - k_d^{(\psi_d)} on common support');
cb = colorbar(ax); cb.Label.String = 'A/A';
ax = nexttile(t);
finite = rows & isfinite(local.L(:, 2, 1));
scatter(ax, x(finite, 1), x(finite, 2), 18, 1000*local.L(finite, 2, 1), 'filled');
hold(ax, 'on');
masked = rows & local.small_Lqd;
scatter(ax, x(masked, 1), x(masked, 2), 40, 'k', 'o');
title(ax, sprintf('L_{qd}; black circles: |L_{qd}| < %.3g mH', 1000*local.denominator_limits(2)));
cb = colorbar(ax); cb.Label.String = 'L_{qd} in mH';
ax = nexttile(t);
finite = rows & local.valid_local;
scatter(ax, x(finite, 1), x(finite, 2), 18, 1000*local.L(finite, 1, 1), 'filled');
title(ax, 'L_{dd}: exploratory saturation association');
cb = colorbar(ax); cb.Label.String = 'L_{dd} in mH';
for ax = findobj(t, 'Type', 'axes')'
    xlabel(ax, 'i_d in A'); ylabel(ax, 'i_q in A'); grid(ax, 'on');
end
ax = nexttile(t);
g = local.component_group_summary;
plot(ax, g(:, 1), g(:, 4:6), '.-', 'LineWidth', 1.2);
hold(ax, 'on');
yline(ax, local.global_contour_k, '--k');
yline(ax, local.global_derivative_k, ':k');
xlabel(ax, 'Nominal i_e in A'); ylabel(ax, 'Factor in A/A'); grid(ax, 'on');
title(ax, 'Group medians: identical support for three local factors');
legend(ax, {'k_d^{(\psi_d)}', 'k_d^{(\psi_q)}', 'k_{LS}', ...
    'Global contour', 'Global derivative'}, 'Location', 'best');
title(t, sprintf('Set%d, nominal i_e = %g A: denominator masks, no clipping or extrapolation', ...
    local.set, local.excitation(group)));
set(findall(plot_tab,'Type','axes'),'Colormap',turbo(256));
if ~isempty(export_format)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'matlab', 'common'));
    export_investigation_plot(t, fullfile(root, 'matlab', 'k_factor_investigation', 'results', ...
        sprintf('local_comparison_set%d_group%02d', local.set, group)), export_format);
end
end
