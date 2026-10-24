%% EESM -> PSM: compare d- and q-axis shifts from constant flux
% kd = -Delta Id / Delta Ie at fixed Iq; kq = -Delta Iq / Delta Ie at fixed Id.
% Endpoints are symmetric about (id_grid, iq_grid). Use the nearest root.
clear; clc;
root = fileparts(mfilename('fullpath'));
data_file = fullfile(root, 'Resources', 'DataSets', 'Set2', 'extracted', ...
    'evaluation', 'eval_flux_eesm_excitation_0A_to_40A_measurements.mat');
speed_field = 'pos_3000_rpm';
n_mesh = 41;
pair_to_plot = 8;  % Figures 1 and 2: excitation pair 20 to 23 in A.
output_dir = fullfile(root, 'results');

%% Load measurements and sort by excitation current
loaded = load(data_file, 'Total');
slices = loaded.Total.(speed_field);
names = fieldnames(slices);
names = names(~cellfun('isempty', regexp(names, '^(pos|neg)_\d+_A$', 'once')));
ie_ref = cellfun(@(name) slices.(name).Ref_exc_current, names);
[ie_ref, order] = sort(ie_ref); names = names(order);
n_exc = numel(names);
current = cell(n_exc, 1); flux = current; triangles = current; edges = current;
for g = 1:n_exc
    meas = slices.(names{g}).Meas;
    x = [meas.PE_Id_meas_avg(:), meas.PE_Iq_meas_avg(:), meas.I_exc_meas(:)];
    psi = [meas.PE_psi_d_meas(:), meas.PE_psi_q_meas(:)];
    keep = all(isfinite([x, psi]), 2);
    x = x(keep, :); psi = psi(keep, :);
    tri = delaunay(x(:, 1), x(:, 2));
    edge = unique(sort([tri(:, [1 2]); tri(:, [2 3]); tri(:, [1 3])], 2), 'rows');
    length_squared = sum((x(edge(:, 1), 1:2)-x(edge(:, 2), 1:2)).^2, 2);
    limit = 9 * mean(length_squared);  % MeasEval long-edge filter.
    edges{g} = edge(length_squared < limit, :);
    triangle_lengths = [sum((x(tri(:, 1), 1:2)-x(tri(:, 2), 1:2)).^2, 2), ...
        sum((x(tri(:, 2), 1:2)-x(tri(:, 3), 1:2)).^2, 2), ...
        sum((x(tri(:, 1), 1:2)-x(tri(:, 3), 1:2)).^2, 2)];
    triangles{g} = tri(all(triangle_lengths < limit, 2), :);
    current{g} = x; flux{g} = psi;
end
all_current = vertcat(current{:}); all_flux = vertcat(flux{:});
id_axis = linspace(min(all_current(:, 1)), max(all_current(:, 1)), n_mesh);
iq_axis = linspace(min(all_current(:, 2)), max(all_current(:, 2)), n_mesh);
[id_grid, iq_grid] = meshgrid(id_axis, iq_axis);

%% Original MeasEval factor, using the middle 25 to 75 percent of psi_d
psi_range = max(all_flux(:, 1))-min(all_flux(:, 1));
legacy_levels = linspace(min(all_flux(:, 1))+.25*psi_range, ...
    min(all_flux(:, 1))+.75*psi_range, 100);
contour_means = nan(n_exc, 100, 3);
for g = 1:n_exc
    x = current{g}; psi = flux{g}; a = edges{g}(:, 1); b = edges{g}(:, 2);
    for p = 1:numel(legacy_levels)
        t = (legacy_levels(p)-psi(a, 1)) ./ (psi(b, 1)-psi(a, 1));
        hit = t >= 0 & t <= 1;
        points = x(a(hit), :) + t(hit).*(x(b(hit), :)-x(a(hit), :));
        contour_means(g, p, :) = mean(points, 1);
    end
end
legacy_ratios = abs(diff(contour_means(:, :, 1), 1, 1)) ./ diff(contour_means(:, :, 3), 1, 1);
legacy_ratios(~isfinite(legacy_ratios) | legacy_ratios < eps) = nan;
k_legacy = mean(legacy_ratios(:), 'omitnan');
fprintf('Original MeasEval factor: %.12f\n', k_legacy);

%% Find equal-flux endpoints in adjacent excitation slices
results = cell(2, 2);  % {shift axis, flux component}
for axis = 1:2
    for component = 1:2
        results{axis, component} = struct('k', nan(n_mesh, n_mesh, n_exc-1), ...
            'psi_const', nan(n_mesh, n_mesh, n_exc-1), ...
            'ie', nan(n_mesh, n_mesh, n_exc-1), ...
            'delta_ie', nan(n_mesh, n_mesh, n_exc-1));
    end
end
for g = 1:n_exc-1
    support_a = triangulation(triangles{g}, current{g}(:, 1:2));
    support_b = triangulation(triangles{g+1}, current{g+1}(:, 1:2));
    for axis = 1:2
        fixed_axis = 3-axis;
        if axis == 1, cuts = iq_axis; centers = id_axis;
        else, cuts = id_axis; centers = iq_axis; end
        for c = 1:n_mesh
            left = line_section(current{g}, flux{g}, edges{g}, fixed_axis, cuts(c));
            right = line_section(current{g+1}, flux{g+1}, edges{g+1}, fixed_axis, cuts(c));
            if size(left, 1) < 2 || size(right, 1) < 2, continue; end
            if axis == 1
                index = sub2ind([n_mesh n_mesh n_exc-1], c*ones(n_mesh, 1), (1:n_mesh)', g*ones(n_mesh, 1));
            else
                index = sub2ind([n_mesh n_mesh n_exc-1], (1:n_mesh)', c*ones(n_mesh, 1), g*ones(n_mesh, 1));
            end
            for component = 1:2
                [k, p, ie, die] = pair_factors(left, right, centers, component+1, ...
                    axis, cuts(c), support_a, support_b);
                results{axis, component}.k(index) = k;
                results{axis, component}.psi_const(index) = p;
                results{axis, component}.ie(index) = ie;
                results{axis, component}.delta_ie(index) = die;
            end
        end
    end
    fprintf('Excitation pair: %g to %g in A\n', ie_ref(g), ie_ref(g+1));
end
kd_psid = results{1, 1}; kd_psiq = results{1, 2};
kq_psid = results{2, 1}; kq_psiq = results{2, 2};
% Array dimensions: q-axis position, d-axis position, excitation pair.
kd_diff = kd_psid.k-kd_psiq.k;
kq_diff = kq_psid.k-kq_psiq.k;
%% Figures 1 and 2: compare at the same mean Id and Iq
sets = {kd_psid, kd_psiq, kq_psid, kq_psiq};
plots = {'kd_comparison.png', 'kq_comparison.png', 'k_over_flux.png', 'flux_original_transformed.png'};
layouts = cell(1, 4);
for f = 1:2
    fig = figure(f); clf(fig); theme(fig, 'light');
    set(fig, 'Name', sprintf('Factor comparison %g', f), 'Position', [40 100 1550 550]);
    layouts{f} = tiledlayout(1, 3, 'TileSpacing', 'compact');
    d = sets{2*f-1}.k(:, :, pair_to_plot); q = sets{2*f}.k(:, :, pair_to_plot);
    values = {d, q, d-q}; factor = {'k_d', 'k_q'};
    titles = {[factor{f} ' from constant \psi_d'], [factor{f} ' from constant \psi_q'], 'Difference: first minus second'};
    for j = 1:3
        nexttile; mesh(id_grid, iq_grid, values{j}, 'FaceColor', 'none');
        title(titles{j}); xlabel('Mean I_d in A'); ylabel('Mean I_q in A');
        zlabel([factor{f} ' in A/A']); colorbar; grid on;
    end
    title(layouts{f}, sprintf('Excitation current: %g to %g in A', ie_ref(pair_to_plot), ie_ref(pair_to_plot+1)));
end

%% Figure 3: all factors as waterfalls, one colour per excitation pair
fig = figure(3); clf(fig); theme(fig, 'light');
set(fig, 'Name', 'Factors over flux', 'Position', [60 60 1550 1000]);
layouts{3} = tiledlayout(2, 2, 'TileSpacing', 'compact');
colors = turbo(n_exc-1);
labels = arrayfun(@(g) sprintf('Ie: %g to %g in A; Delta Ie: %g in A', ...
    ie_ref(g), ie_ref(g+1), ie_ref(g+1)-ie_ref(g)), 1:n_exc-1, 'UniformOutput', false);
titles = {'k_d from constant \psi_d', 'k_d from constant \psi_q', ...
    'k_q from constant \psi_d', 'k_q from constant \psi_q'};
for j = 1:4
    nexttile; hold on;
    for g = 1:n_exc-1
        p = sets{j}.psi_const(:, :, g); k = sets{j}.k(:, :, g); ie = sets{j}.ie(:, :, g);
        if j > 2, p = p.'; k = k.'; ie = ie.'; end
        waterfall(p, ie, k, 'EdgeColor', colors(g, :), 'FaceColor', 'none', 'LineWidth', .5);
    end
    title(titles{j}); xlabel('Constant flux in Wb'); ylabel('Mean excitation current in A');
    zlabel('k in A/A'); view(40, 25); grid on;
    legend(labels, 'Location', 'eastoutside', 'FontSize', 7, 'Interpreter', 'none');
end

%% Figure 4: colour identifies excitation current, before and after transformation
fig = figure(4); clf(fig); theme(fig, 'light');
set(fig, 'Name', 'Original and transformed flux maps', 'Position', [80 80 1500 950]);
layouts{4} = tiledlayout(2, 2, 'TileSpacing', 'compact');
colors = turbo(n_exc);
labels = arrayfun(@(ie) sprintf('Ie = %g in A', ie), ie_ref, 'UniformOutput', false);
flux_labels = {'\psi_d in Wb', '\psi_q in Wb'};
for component = 1:2
    for transformed = 0:1
        nexttile; hold on;
        for g = 1:n_exc
            x = current{g}; d = x(:, 1)+transformed*k_legacy*x(:, 3);
            trisurf(triangles{g}, d, x(:, 2), flux{g}(:, component), ...
                'FaceColor', 'none', 'EdgeColor', colors(g, :), 'LineWidth', .5);
        end
        if transformed
            xlabel('I_d^* = I_d + k I_e in A'); title(sprintf('Transformed, k = %.5f', k_legacy));
        else
            xlabel('d-axis current in A'); title('Original');
        end
        ylabel('q-axis current in A'); zlabel(flux_labels{component}); view(3); grid on;
        legend(labels, 'Location', 'eastoutside', 'FontSize', 7, 'Interpreter', 'none');
    end
end

%% Save the calculated factors and figures
if ~isfolder(output_dir), mkdir(output_dir); end
save(fullfile(output_dir, 'winding_factors.mat'), 'ie_ref', 'id_grid', 'iq_grid', ...
    'kd_psid', 'kd_psiq', 'kq_psid', 'kq_psiq', 'kd_diff', 'kq_diff', 'k_legacy');
for f = 1:4
    exportgraphics(layouts{f}, fullfile(output_dir, plots{f}), 'Resolution', 150);
end
fprintf('MAT file and figures saved to %s\n', output_dir);

%% Two numerical helpers
function section = line_section(x, psi, edges, fixed_axis, cut)
% Intersect Delaunay edges with a constant Id or Iq line.
    a = edges(:, 1); b = edges(:, 2); values = [x, psi];
    t = (cut-x(a, fixed_axis)) ./ (x(b, fixed_axis)-x(a, fixed_axis));
    hit = t >= 0 & t <= 1;
    points = values(a(hit), :) + t(hit).*(values(b(hit), :)-values(a(hit), :));
    points = [points; values(abs(x(:, fixed_axis)-cut) < 1e-9, :)];
    [~, keep] = unique(round(points(:, 3-fixed_axis), 8));
    section = points(keep, [3-fixed_axis, 4, 5, 3]); % free current, psi_d, psi_q, Ie
end

function [k, psi, ie, die] = pair_factors(a, b, centers, component, axis, cut, support_a, support_b)
    k = nan(numel(centers), 1); psi = k; ie = k; die = k;
    for p = 1:numel(centers)
        center = centers(p);
        low = max(2*(a(1, 1)-center), 2*(center-b(end, 1)));
        high = min(2*(a(end, 1)-center), 2*(center-b(1, 1)));
        if low >= high, continue; end
        knots = unique([low; high; 2*(a(:, 1)-center); 2*(center-b(:, 1))]);
        knots = knots(knots >= low & knots <= high);
        difference = interp1(a(:, 1), a(:, component), center+knots/2) ...
            - interp1(b(:, 1), b(:, component), center-knots/2);
        cross = find(difference(1:end-1).*difference(2:end) <= 0 ...
            & difference(1:end-1) ~= difference(2:end));
        delta = knots(cross)-difference(cross).*(knots(cross+1)-knots(cross)) ...
            ./ (difference(cross+1)-difference(cross));
        if isempty(delta), continue; end
        points_a = zeros(numel(delta), 2); points_b = points_a;
        points_a(:, axis) = center+delta/2; points_b(:, axis) = center-delta/2;
        points_a(:, 3-axis) = cut; points_b(:, 3-axis) = cut;
        inside = ~isnan(pointLocation(support_a, points_a)) & ~isnan(pointLocation(support_b, points_b));
        delta = delta(inside);
        if isempty(delta), continue; end
        [~, nearest] = min(abs(delta)); delta = delta(nearest);
        ie_a = interp1(a(:, 1), a(:, 4), center+delta/2);
        ie_b = interp1(b(:, 1), b(:, 4), center-delta/2);
        if abs(ie_b-ie_a) < 1e-3, continue; end
        k(p) = delta/(ie_b-ie_a);
        psi(p) = interp1(a(:, 1), a(:, component), center+delta/2);
        ie(p) = (ie_a+ie_b)/2; die(p) = ie_b-ie_a;
    end
end

