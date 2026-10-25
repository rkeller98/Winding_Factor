function [reduction, fig] = plot_contour_factor(currents, flux, groups, k, contour, excitation)
% Display the contour estimate and its field reduction in one tabbed figure.
if nargin < 6
    excitation = [];
end

matlab_folder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(matlab_folder, 'common'));

current_scale = max(abs(currents), [], 1);
map = fit_reduction_map(currents, flux, 1, k, current_scale, 81, true);
reconstructed_flux = reduction_lookup(map, currents);
flux_error = reconstructed_flux - flux;
error_limit = max(1000 * abs(flux_error), [], 'all', 'omitnan');
error_limit = max(error_limit, eps);

reduction.s = currents(:, 1) + k * currents(:, 3);
reduction.map = map;
reduction.psi_reconstructed = reconstructed_flux;
reduction.error = flux_error;

group_ids = unique(groups);
triangles = cell(numel(group_ids), 1);

for group_index = 1:numel(group_ids)
    group_rows = groups == group_ids(group_index);
    triangles{group_index} = measurement_triangles( ...
        currents(group_rows, 1:2), current_scale(1:2));
end
reduction.triangles = triangles;

fig = figure('Visible', 'on', 'Name', 'Contour-factor analysis', ...
    'Tag', 'contour_factor', 'Position', [30 30 1500 900]);
theme(fig, 'light');
tabs = uitabgroup(fig, 'Tag', 'contour_tabs');

measured_tab = uitab(tabs, 'Title', 'Measured flux maps', 'Tag', 'measured_flux');
shifted_tab = uitab(tabs, 'Title', 'Shifted measured maps', 'Tag', 'shifted_flux');
reduced_tab = uitab(tabs, 'Title', 'Fixed reduced reference map', 'Tag', 'reduced_flux');
reconstructed_tab = uitab(tabs, 'Title', 'Reconstruction', 'Tag', 'reconstructed_flux');
comparison_tab = uitab(tabs, 'Title', 'Measured + reconstructed', 'Tag', 'measured_reconstructed_flux');
error_tab = uitab(tabs, 'Title', 'Reconstruction error');
error_tabs = uitabgroup(error_tab);
excitation_limits = [min(currents(:, 3)), max(currents(:, 3))];
excitation_fraction = (currents(:, 3) - excitation_limits(1)) ...
    / max(diff(excitation_limits), eps);
measured_colors = interp1(linspace(0, 1, 256), turbo(256), excitation_fraction);

% Keep each group's original mesh, also when shifting its d-current to s.
field_values = {flux, flux, reconstructed_flux, flux};
field_names = {'Measured flux fields', ...
    sprintf('Measured fields shifted to s = i_d + k i_e, k = %.5f A/A', k), ...
    'Flux reconstructed from F(s, i_q)', 'Measured and reconstructed flux fields'};
field_tabs = [measured_tab, shifted_tab, reconstructed_tab, comparison_tab];
current_columns = {[2 1], [1 2]};
current_labels = {'i_q in A', 'i_d in A'; 'i_d in A', 'i_q in A'};
field_meshes = gobjects(numel(group_ids), 2, numel(field_tabs));
comparison_meshes = gobjects(numel(group_ids), 2);
comparison_points = gobjects(numel(group_ids), 2);
group_checks = gobjects(numel(group_ids), numel(field_tabs));
show_buttons = gobjects(numel(field_tabs), 1);
hide_buttons = gobjects(numel(field_tabs), 1);

for field_index = 1:numel(field_tabs)
    plot_tab = field_tabs(field_index);
    plot_panel = uipanel(plot_tab, 'BorderType', 'none', ...
        'Units', 'normalized', 'Position', [0 0 0.86 1]);
    layout = tiledlayout(plot_panel, 1, 2, 'TileSpacing', 'compact');
    selection_panel = uipanel(plot_tab, 'Title', 'Excitation groups', ...
        'FontSize', 11, 'Units', 'normalized', 'Position', [0.86 0.02 0.14 0.96]);
    show_buttons(field_index) = uicontrol(selection_panel, 'Style', 'pushbutton', ...
        'String', 'All', 'Tag', 'show_all_fields', ...
        'Units', 'normalized', 'Position', [0.05 0.92 0.42 0.05]);
    hide_buttons(field_index) = uicontrol(selection_panel, 'Style', 'pushbutton', ...
        'String', 'None', 'Tag', 'hide_all_fields', ...
        'Units', 'normalized', 'Position', [0.53 0.92 0.42 0.05]);
    row_height = min(0.05, 0.75 / numel(group_ids));
    for group_index = 1:numel(group_ids)
        group_name = sprintf('Group %d', group_index);
        if ~isempty(excitation)
            group_name = sprintf('%g A (reference)', excitation(group_index));
        end
        group_checks(group_index, field_index) = uicontrol(selection_panel, ...
            'Style', 'checkbox', 'String', group_name, 'Value', 1, 'FontSize', 11, ...
            'Tag', sprintf('field_group%02d', group_index), ...
            'Units', 'normalized', ...
            'Position', [0.05 0.87 - group_index * row_height 0.9 row_height]);
    end
    uicontrol(selection_panel, 'Style', 'text', 'FontSize', 10, ...
        'String', {'Selection applies to', 'all flux-map tabs.'}, ...
        'Units', 'normalized', 'Position', [0.05 0.01 0.9 0.08]);
    values = field_values{field_index};
    plot_currents = currents;
    axis_labels = current_labels;

    if field_index == 2
        plot_currents(:, 1) = reduction.s;
        axis_labels{1, 2} = 's = i_d + k i_e in A';
        axis_labels{2, 1} = 's = i_d + k i_e in A';
    end

    for component = 1:2
        nexttile(layout);
        hold on;

        for group_index = 1:numel(group_ids)
            group_rows = groups == group_ids(group_index);
            group_currents = plot_currents(group_rows, :);
            group_flux = values(group_rows, component);
            mesh = triangles{group_index};
            mesh = mesh(all(isfinite(group_flux(mesh)), 2), :);

            columns = current_columns{component};
            field_meshes(group_index, component, field_index) = trisurf( ...
                mesh, group_currents(:, columns(1)), group_currents(:, columns(2)), ...
                1000 * group_flux, group_currents(:, 3), ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 0.5);

            if plot_tab == comparison_tab
                set(field_meshes(group_index, component, field_index), ...
                    'FaceColor', 'interp', 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
                    'FaceVertexCData', measured_colors(group_rows, :));
                reconstructed_values = reconstructed_flux(group_rows, component);
                reconstructed_mesh = mesh(all(isfinite(reconstructed_values(mesh)), 2), :);
                comparison_meshes(group_index, component) = trisurf( ...
                    reconstructed_mesh, group_currents(:, columns(1)), ...
                    group_currents(:, columns(2)), 1000 * reconstructed_values, ...
                    'FaceColor', 'none', 'EdgeColor', [0.12 0.12 0.12], ...
                    'LineStyle', '-', 'LineWidth', 0.35, 'EdgeAlpha', 0.35);
                point_error = 1000 * flux_error(group_rows, component);
                supported_points = isfinite(reconstructed_values) & isfinite(point_error);
                comparison_points(group_index, component) = scatter3( ...
                    group_currents(supported_points, columns(1)), ...
                    group_currents(supported_points, columns(2)), ...
                    1000 * reconstructed_values(supported_points), ...
                    24, point_error(supported_points), 'filled', 'MarkerEdgeColor', 'none');
            end
        end

        if plot_tab == comparison_tab
            measured_key = plot3(NaN, NaN, NaN, 's', ...
                'MarkerFaceColor', [0.2 0.5 0.9], 'MarkerEdgeColor', 'none', 'MarkerSize', 9);
            reconstructed_key = plot3(NaN, NaN, NaN, '-o', ...
                'Color', [0.12 0.12 0.12], 'LineWidth', 0.6, ...
                'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'none');
            legend([measured_key, reconstructed_key], ...
                {'Measured (surface)', 'Reconstructed (error color)'}, 'Location', 'southoutside', ...
                'Orientation', 'horizontal', 'FontSize', 10, 'AutoUpdate', 'off');
        end

        xlabel(axis_labels{component, 1});
        ylabel(axis_labels{component, 2});
        zlabel(flux_label(component));
        title(flux_label(component));
        view(35, 25);
        grid on;
        axis manual;
        color_scale = colorbar;
        if plot_tab == comparison_tab
            clim([-error_limit error_limit]);
            color_scale.Label.String = 'Reconstructed - measured in mVs';
        else
            clim(excitation_limits);
            color_scale.Label.String = 'Measured i_e in A';
        end
    end

    if plot_tab == comparison_tab
        set(findall(plot_tab, 'Type', 'axes'), 'Colormap', error_colors);
    else
        set(findall(plot_tab, 'Type', 'axes'), 'Colormap', turbo(256));
    end
    title(layout, field_names{field_index});
end

% Include both overlaid fields in the same group selection.
field_meshes = cat(3, field_meshes, comparison_meshes, comparison_points);
% Keep the same group selection across both components and all field tabs.
for field_index = 1:numel(field_tabs)
    show_buttons(field_index).Callback = ...
        @(~, ~) set_field_visibility(field_meshes, group_checks, true);
    hide_buttons(field_index).Callback = ...
        @(~, ~) set_field_visibility(field_meshes, group_checks, false);
    for group_index = 1:numel(group_ids)
        group_checks(group_index, field_index).Callback = @(source, ~) ...
            set_field_visibility(field_meshes(group_index, :, :), ...
            group_checks(group_index, :), source.Value);
    end
end

% The stored reduced field is the same 81-by-81 LUT used by method A.
s_axis = map.origin(1) + map.scale(1) * map.axes{1};
iq_axis = map.origin(2) + map.scale(2) * map.axes{2};
[s_grid, iq_grid] = ndgrid(s_axis, iq_axis);

layout = tiledlayout(reduced_tab, 1, 2, 'TileSpacing', 'compact');

for component = 1:2
    nexttile(layout);
    if component == 1
        surf(iq_grid, s_grid, 1000 * map.values{component}, ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 0.5);
        xlabel('i_q in A');
        ylabel('s = i_d + k i_e in A');
    else
        surf(s_grid, iq_grid, 1000 * map.values{component}, ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 0.5);
        xlabel('s = i_d + k i_e in A');
        ylabel('i_q in A');
    end
    zlabel(flux_label(component));
    title(flux_label(component));
    view(35, 25);
    grid on;
    color_scale = colorbar;
    color_scale.Label.String = flux_label(component);
end

set(findall(reduced_tab, 'Type', 'axes'), 'Colormap', parula(256));
title(layout, sprintf('Fixed reduced reference map F_A(s, i_q); k_A = %.5f A/A; no extrapolation', k));

% Plot signed reconstruction error on each original excitation group's mesh.
% The heatmap is a top view of that mesh, with no extra fitted error field.
for group_index = 1:numel(group_ids)
    group_rows = groups == group_ids(group_index);
    group_currents = currents(group_rows, :);
    group_error = 1000 * flux_error(group_rows, :);
    mesh = triangles{group_index};

    group_title = sprintf('Group %d', group_index);
    if ~isempty(excitation)
        group_title = sprintf('%s: %g A (reference)', group_title, excitation(group_index));
    end
    plot_tab = uitab(error_tabs, ...
        'Title', group_title, ...
        'Tag', sprintf('reconstruction_error_group%02d', group_index));
    layout = tiledlayout(plot_tab, 2, 3, 'TileSpacing', 'compact');

    for component = 1:2
        values = group_error(:, component);
        supported_mesh = mesh(all(isfinite(values(mesh)), 2), :);

        nexttile(layout);
        columns = current_columns{component};
        trisurf(supported_mesh, group_currents(:, columns(1)), group_currents(:, columns(2)), ...
            values, values, 'FaceColor', 'none', 'EdgeColor', 'interp');
        xlabel(current_labels{component, 1});
        ylabel(current_labels{component, 2});
        zlabel('Error in mVs');
        title(sprintf('%s: error mesh', flux_label(component)));
        clim([-error_limit error_limit]);
        colorbar;
        grid on;
        view(35, 25);

        nexttile(layout);
        trisurf(supported_mesh, group_currents(:, 1), group_currents(:, 2), ...
            zeros(size(values)), values, 'FaceColor', 'interp', 'EdgeColor', 'none');
        xlabel('i_d in A');
        ylabel('i_q in A');
        title(sprintf('%s: error heatmap', flux_label(component)));
        clim([-error_limit error_limit]);
        color_scale = colorbar;
        color_scale.Label.String = 'Error in mVs';
        axis tight;
        view(2);

        nexttile(layout);
        error_samples = values(isfinite(values));
        error_mean = mean(error_samples);
        error_std = std(error_samples);
        error_rmse = sqrt(mean(error_samples.^2));

        if ~isempty(error_samples)
            histogram(error_samples, 'Normalization', 'pdf', ...
                'FaceColor', [0.2 0.45 0.75], 'EdgeColor', 'white');
            hold on;
            xline(error_mean, 'k-', '\mu', 'LineWidth', 1.5);
            xline(error_mean - error_std, 'k--', '\mu-\sigma');
            xline(error_mean + error_std, 'k--', '\mu+\sigma');
        else
            text(0.5, 0.5, 'No supported points', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
        end

        xlabel('Error in mVs');
        ylabel('Density in 1/mVs');
        title({sprintf('%s: error distribution', flux_label(component)), ...
            sprintf('\\mu = %.3f, \\sigma = %.3f mVs', error_mean, error_std), ...
            sprintf('RMSE = %.3f mVs | n = %d', error_rmse, numel(error_samples))});
        grid on;
    end

    set(findall(plot_tab, 'Type', 'axes'), 'Colormap', error_colors);
    title(layout, sprintf('Reconstructed - measured | %s', group_title));
end

contour_tab = uitab(tabs, 'Title', 'Contour factor', 'Tag', 'contours');
plot_contour_curves(contour_tab, k, contour, excitation);

q_tab = uitab(tabs, 'Title', '3D diagnostics', 'Tag', 'contour_q_dependence');
plot_q_dependence(q_tab, contour);

tabs.SelectedTab = measured_tab;
end


function set_field_visibility(meshes, checkboxes, show)
visibility = 'off';
if show
    visibility = 'on';
end
set(meshes(:), 'Visible', visibility);
set(checkboxes(:), 'Value', show);
end


function plot_contour_curves(plot_tab, k, contour, excitation)
full = contour.full_range;
flux_levels = 1000 * full.levels;
fit_limits = 1000 * contour.levels([1 end]);
if isempty(excitation)
    group_names = compose('Group %d', (1:size(full.means, 1))');
else
    group_names = compose('%g A', excitation(:));
end
pair_names = group_names(1:end-1) + " \rightarrow " + group_names(2:end);

layout = tiledlayout(plot_tab, 2, 2, 'TileSpacing', 'compact');

values = {
    full.ratios
    full.signed
    full.iq_drift
    1000 * full.centroid_error
};
labels = {
    'Absolute-numerator k in A/A'
    'Signed -\Delta i_d / \Delta i_e in A/A'
    '\Delta i_q / \Delta i_e in A/A'
    'Centroid flux minus contour level in mVs'
};
panel_titles = {
    sprintf('Current estimate: k = %.5f A/A', k)
    'Signed contour secants'
    'q-current drift between contour centroids'
    'Contour-centroid flux departure'
};

for panel = 1:4
    ax = nexttile(layout);
    curves = plot(flux_levels, values{panel}', 'LineWidth', 0.9);
    curve_colors = turbo(numel(curves));
    for curve_index = 1:numel(curves)
        curves(curve_index).Color = curve_colors(curve_index, :);
    end
    hold on;

    if panel == 1
        yline(k, 'k-', sprintf('k = %.5f', k), 'LineWidth', 1.5);
    end

    y_limits = ylim;
    patch(ax, fit_limits([1 2 2 1]), y_limits([1 1 2 2]), [0.4 0.7 1], ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none');
    xline(fit_limits(1), 'k--', '25%');
    xline(fit_limits(2), 'k--', '75%');
    ylim(y_limits);
    xlim(flux_levels([1 end]));
    xlabel('Contour level \psi_d in mVs');
    ylabel(labels{panel});
    title(panel_titles{panel});
    grid on;

    if panel == 4
        curve_names = group_names;
    else
        curve_names = pair_names;
    end
    curve_legend = legend(ax, curves, curve_names, ...
        'Location', 'eastoutside', 'FontSize', 9, 'AutoUpdate', 'off');
    if ~isempty(excitation)
        curve_legend.Title.String = 'i_{e,req}';
    end
end

title(layout, 'Full d-flux value range; only the shaded 25--75% interval determines k');
end


function plot_q_dependence(plot_tab, contour)
% Each point is one existing adjacent-group secant, not a fixed-i_q estimate.
full = contour.full_range;
paired_iq = (full.means(1:end-1, :, 2) + full.means(2:end, :, 2)) / 2;
paired_ie = (full.means(1:end-1, :, 3) + full.means(2:end, :, 3)) / 2;
flux_levels = 1000 * full.levels;

layout = tiledlayout(plot_tab, 1, 2, 'TileSpacing', 'compact');
values = {full.ratios, full.iq_drift};
labels = {'Contour k in A/A', '\Delta i_q / \Delta i_e in A/A'};
excitation_limits = [min(paired_ie, [], 'all', 'omitnan'), ...
    max(paired_ie, [], 'all', 'omitnan')];

for panel = 1:2
    nexttile(layout);
    hold on;

    for pair_index = 1:size(paired_iq, 1)
        iq = paired_iq(pair_index, :);
        ie = paired_ie(pair_index, :);
        value = values{panel}(pair_index, :);
        supported = isfinite(iq) & isfinite(ie) & isfinite(value);

        scatter3(flux_levels(supported), iq(supported), value(supported), ...
            18, ie(supported), 'filled');
    end

    xlabel('\psi_d in mVs');
    ylabel('Mean contour-centroid i_q in A');
    zlabel(labels{panel});
    title(labels{panel});
    grid on;
    view(40, 30);
    clim(excitation_limits);
    color_scale = colorbar;
    color_scale.Label.String = 'Mean i_e in A';
end

set(findall(plot_tab, 'Type', 'axes'), 'Colormap', turbo(256));
title(layout, 'Contour-centroid diagnostics: k and q-current drift after contour averaging');
end


function triangles = measurement_triangles(id_iq, current_scale)
% Keep the existing original-surface mesh and long-edge cutoff.
coordinates = id_iq ./ current_scale;
mesh = delaunayTriangulation(coordinates);
triangles = mesh.ConnectivityList;
edges = mesh.edges;

edge_length_squared = sum( ...
    (coordinates(edges(:, 1), :) - coordinates(edges(:, 2), :)).^2, 2);
length_limit = 9 * mean(edge_length_squared);
supported = true(size(triangles, 1), 1);
vertex_pairs = [1 2; 1 3; 2 3];

for edge_index = 1:size(vertex_pairs, 1)
    pair = vertex_pairs(edge_index, :);
    edge_vectors = coordinates(triangles(:, pair(1)), :) ...
        - coordinates(triangles(:, pair(2)), :);
    supported = supported & sum(edge_vectors.^2, 2) < length_limit;
end

triangles = triangles(supported, :);
end


function label = flux_label(component)
labels = {'\psi_d in mVs', '\psi_q in mVs'};
label = labels{component};
end


function colors = error_colors
colors = [
    linspace(0.1, 1, 128)', linspace(0.3, 1, 128)', ones(128, 1)
    ones(128, 1), linspace(1, 0.2, 128)', linspace(1, 0.1, 128)'
];
end
