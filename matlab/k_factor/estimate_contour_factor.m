function [k, detail, fig] = estimate_contour_factor(currents, flux, groups, make_plots, excitation)
% Estimate k from contour intersections in adjacent excitation groups.
% currents = [id_meas, iq_meas, ie_meas], flux = [psi_d, psi_q].
% Without output arguments, show the fields and diagnostics in one tabbed figure.
% Optional excitation contains nominal group currents, in sorted group order.

if nargin < 4
    make_plots = nargout == 0;
end
if nargin < 5
    excitation = [];
end
fig = gobjects(0);

% Use 100 contour levels in the middle 25--75% of the pooled d-flux range.
d_flux = flux(:, 1);
flux_span = max(d_flux) - min(d_flux);
levels = linspace( ...
    min(d_flux) + 0.25 * flux_span, ...
    min(d_flux) + 0.75 * flux_span, 100);
fit_level_count = numel(levels);

% Additional levels are for display only and never enter the estimated mean.
if make_plots
    full_range_levels = linspace(min(d_flux), max(d_flux), 201);
    levels = [levels, full_range_levels];
end

group_ids = unique(groups);
group_count = numel(group_ids);
contour_means = nan(group_count, numel(levels), 3);
centroid_error = nan(group_count, numel(levels));

for group_index = 1:group_count
    group_rows = groups == group_ids(group_index);
    group_currents = currents(group_rows, :);
    group_flux = flux(group_rows, :);

    % Build each excitation group's mesh in measured id/iq coordinates.
    triangles = delaunay(group_currents(:, 1), group_currents(:, 2));
    triangle_edges = [
        triangles(:, [1 2])
        triangles(:, [1 3])
        triangles(:, [2 3])
    ];
    edges = unique(sort(triangle_edges, 2), 'rows');

    edge_vectors = group_currents(edges(:, 1), 1:2) ...
        - group_currents(edges(:, 2), 1:2);
    edge_length_squared = sum(edge_vectors.^2, 2);
    short_edges = edge_length_squared < 9 * mean(edge_length_squared);
    edges = edges(short_edges, :);

    start_rows = edges(:, 1);
    end_rows = edges(:, 2);
    d_flux_at_centroid = scatteredInterpolant( ...
        group_currents(:, 1:2), group_flux(:, 1), 'linear', 'none');

    % Interpolate all three measured currents at each contour crossing.
    for level_index = 1:numel(levels)
        weights = (levels(level_index) - group_flux(start_rows, 1)) ...
            ./ (group_flux(end_rows, 1) - group_flux(start_rows, 1));
        crossing_edges = weights >= 0 & weights <= 1;

        start_currents = group_currents(start_rows(crossing_edges), :);
        end_currents = group_currents(end_rows(crossing_edges), :);
        intersections = start_currents ...
            + weights(crossing_edges) .* (end_currents - start_currents);

        contour_means(group_index, level_index, :) = mean(intersections, 1);

        centroid_currents = reshape( ...
            contour_means(group_index, level_index, 1:2), 1, 2);
        centroid_error(group_index, level_index) = ...
            d_flux_at_centroid(centroid_currents) - levels(level_index);
    end
end

% Adjacent groups share the same selected d-flux level.
delta_id = diff(contour_means(:, :, 1));
delta_ie = diff(contour_means(:, :, 3));
signed_ratios = -delta_id ./ delta_ie;
ratios = abs(delta_id) ./ delta_ie;

invalid_ratios = ~isfinite(ratios) | abs(ratios) < eps;
ratios(invalid_ratios) = NaN;
k = mean(ratios(:, 1:fit_level_count), 'all', 'omitnan');

iq_drift = diff(contour_means(:, :, 2)) ./ delta_ie;
fit_columns = 1:fit_level_count;
detail.levels = levels(fit_columns);
detail.means = contour_means(:, fit_columns, :);
detail.ratios = ratios(:, fit_columns);
detail.signed = signed_ratios(:, fit_columns);
detail.centroid_error = centroid_error(:, fit_columns);
detail.iq_drift = iq_drift(:, fit_columns);

if make_plots
    plot_columns = fit_level_count + (1:numel(full_range_levels));
    detail.full_range.levels = levels(plot_columns);
    detail.full_range.means = contour_means(:, plot_columns, :);
    detail.full_range.ratios = ratios(:, plot_columns);
    detail.full_range.signed = signed_ratios(:, plot_columns);
    detail.full_range.centroid_error = centroid_error(:, plot_columns);
    detail.full_range.iq_drift = iq_drift(:, plot_columns);

    [detail.reduction, fig] = plot_contour_factor( ...
        currents, flux, groups, k, detail, excitation);
end
end
