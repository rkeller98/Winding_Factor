function figures = k_factor_contour_analysis_plot(meas, meta, selected_idx)
% Visualize flux maps, Delaunay contours, centroids and winding factors.

if nargin < 3 || isempty(selected_idx)
    selected_idx = round(linspace(0.1, 0.9, 5) * ...
        (numel(meta.psi_d_values) - 1)) + 1;
end

selected_idx = unique(selected_idx, 'stable');
colors = lines(numel(selected_idx));

N_ie = numel(meta.ie_op);

%% Figure and tabs

figures = figure( ...
    'Name', 'K-Factor Analysis', ...
    'Position', [50 50 1600 800]);

tabs = uitabgroup(figures);

tab_maps = uitab(tabs, ...
    'Title', 'Flux Maps');

tab_error = uitab(tabs, ...
    'Title', 'Centroid Error');

tab_kde = uitab(tabs, ...
    'Title', 'k_de');

ie_tabs = uitabgroup(tab_maps);

% Common color limits across excitation currents
pd_all = meas.psi_d(isfinite(meas.psi_d));
pq_all = meas.psi_q(isfinite(meas.psi_q));

pd_limits = [min(pd_all), max(pd_all)];
pq_limits = [min(pq_all), max(pq_all)];

%% Flux maps and contours

for idx_ie = 1:N_ie

    % Exactly the points used in the calculation
    rows = meta.valid_rows{idx_ie};

    id = meas.id_meas(rows); id = id(:);
    iq = meas.iq_meas(rows); iq = iq(:);
    pd = meas.psi_d(rows); pd = pd(:);
    pq = meas.psi_q(rows); pq = pq(:);

    DT = meta.triangles{idx_ie};

    tab = uitab(ie_tabs, ...
        'Title', sprintf('%g A', meta.ie_op(idx_ie)));

    layout = tiledlayout(tab, 1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %% Psi_d in 2D

    ax1 = nexttile(layout);

    trisurf(DT, id, iq, zeros(size(id)), pd, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');

    view(ax1, 2);
    hold(ax1, 'on');
    axis(ax1, 'equal');
    grid(ax1, 'on');
    colormap(ax1, turbo(256));

    if diff(pd_limits) > 0
        clim(ax1, pd_limits);
    end

    colorbar(ax1);

    xlabel(ax1, 'i_d [A]');
    ylabel(ax1, 'i_q [A]');
    title(ax1, '\psi_d: contours and centroids');

    %% Psi_q in 2D

    ax2 = nexttile(layout);

    DT_q = DT(all(isfinite(pq(DT)), 2), :);

    if ~isempty(DT_q)
        trisurf(DT_q, id, iq, zeros(size(id)), pq, ...
            'EdgeColor', 'none', ...
            'FaceColor', 'interp');

        colorbar(ax2);
    end

    view(ax2, 2);
    axis(ax2, 'equal');
    grid(ax2, 'on');
    colormap(ax2, turbo(256));

    if diff(pq_limits) > 0
        clim(ax2, pq_limits);
    end

    xlabel(ax2, 'i_d [A]');
    ylabel(ax2, 'i_q [A]');
    title(ax2, '\psi_q');

    %% Psi_d in 3D

    ax3 = nexttile(layout);

    trisurf(DT, id, iq, pd, pd, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp', ...
        'FaceAlpha', 0.65);

    hold(ax3, 'on');
    grid(ax3, 'on');
    view(ax3, 35, 25);
    colormap(ax3, turbo(256));

    if diff(pd_limits) > 0
        clim(ax3, pd_limits);
    end

    colorbar(ax3);

    xlabel(ax3, 'i_d [A]');
    ylabel(ax3, 'i_q [A]');
    zlabel(ax3, '\psi_d [Vs]');
    title(ax3, '\psi_d: contour vs. centroid flux');

    %% Selected contour levels

    for idx = 1:numel(selected_idx)

        idx_psi = selected_idx(idx);
        psi_level = meta.psi_d_values(idx_psi);
        points = meta.contour_points{idx_ie, idx_psi};

        if isempty(points)
            continue
        end

        color = colors(idx, :);

        id_c = meta.mean_id(idx_ie, idx_psi);
        iq_c = meta.mean_iq(idx_ie, idx_psi);
        psi_c = meta.psi_at_centroid(idx_ie, idx_psi);

        % Delaunay intersections in 2D
        scatter(ax1, points(:,1), points(:,2), ...
            12, color, 'filled');

        % Delaunay intersections in 3D
        scatter3(ax3, points(:,1), points(:,2), ...
            points(:,3), 12, color, 'filled');

        % Centroid in 2D
        plot(ax1, id_c, iq_c, 'p', ...
            'MarkerSize', 12, ...
            'MarkerFaceColor', color, ...
            'MarkerEdgeColor', 'k');

        % Centroid at nominal contour level
        plot3(ax3, id_c, iq_c, psi_level, 'o', ...
            'MarkerSize', 9, ...
            'MarkerFaceColor', color, ...
            'MarkerEdgeColor', 'k');

        if isfinite(psi_c)

            % Interpolated flux at centroid
            plot3(ax3, id_c, iq_c, psi_c, 'kx', ...
                'MarkerSize', 10, ...
                'LineWidth', 2);

            % Difference between contour and centroid flux
            plot3(ax3, ...
                [id_c id_c], ...
                [iq_c iq_c], ...
                [psi_level psi_c], '-', ...
                'Color', color, ...
                'LineWidth', 2);

            % Numerical centroid flux in mVs
            text(ax1, id_c, iq_c, ...
                sprintf(' %.1f', 1000 * psi_c), ...
                'Color', color, ...
                'FontSize', 9);

        end
    end
end

%% Centroid flux deviation

layout_error = tiledlayout(tab_error, 1, 1, ...
    'Padding', 'compact');

ax = nexttile(layout_error);

plot(ax, ...
    1000 * meta.psi_d_values, ...
    1000 * meta.centroid_error.', ...
    'LineWidth', 1.2);

xlabel(ax, '\psi_d [mVs]');
ylabel(ax, '\Delta\psi_d [mVs]');
title(ax, 'Flux deviation at contour centroids');

legend(ax, compose('%g A', meta.ie_op), ...
    'Location', 'eastoutside');

grid(ax, 'on');

%% k_de contour analysis

layout_kde = tiledlayout(tab_kde, 1, 1, ...
    'Padding', 'compact');

ax = nexttile(layout_kde);
hold(ax, 'on');

for idx = 1:size(meta.ie_pairs, 1)

    plot(ax, ...
        1000 * meta.psi_d_values, ...
        meta.k_signed(idx,:), ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('%g A -> %g A', ...
        meta.ie_pairs(idx,1), ...
        meta.ie_pairs(idx,2)));

end

hold(ax, 'off');

xlabel(ax, '\psi_d [mVs]');
ylabel(ax, 'k_{de} [A/A]');
title(ax, 'Contour-based winding factor');

legend(ax, 'show', 'Location', 'eastoutside');

grid(ax, 'on');

%% Common plot appearance

all_axes = findall(figures, 'Type', 'axes');

set(all_axes, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1, ...
    'Box', 'on', ...
    'GridAlpha', 0.15);

tabs.SelectedTab = tab_maps;

end