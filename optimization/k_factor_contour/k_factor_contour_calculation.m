function [kde, meta] = k_factor_contour_calculation(meas)
% Conventional Delaunay contour factor and data for subsequent analysis.
% One speed per dataset. Currents in A, flux linkages in Vs.
% kde matches the legacy 25-75% contour estimator (100 contour levels).

id = meas.id_meas(:);
iq = meas.iq_meas(:);
ie = meas.ie_meas(:);
ie_req = meas.ie_req(:);
psi_d = meas.psi_d(:);

ie_op = unique(ie_req(isfinite(ie_req)));
psi_valid = psi_d(isfinite(psi_d));
psi_min = min(psi_valid);
psi_span = max(psi_valid) - psi_min;

% Full range for diagnostics, plus exactly the legacy levels for kde.
psi_full = linspace(psi_min, psi_min + psi_span, 201);
psi_legacy = linspace(psi_min + 0.25 * psi_span, ...
                      psi_min + 0.75 * psi_span, 100);
psi_d_values = unique([psi_full, psi_legacy]);
legacy_mask = ismember(psi_d_values, psi_legacy);

N_ie = numel(ie_op);
N_psi = numel(psi_d_values);

contour_points = cell(N_ie, N_psi);
mean_id = nan(N_ie, N_psi);
mean_iq = nan(N_ie, N_psi);
mean_ie = nan(N_ie, N_psi);
num_points = zeros(N_ie, N_psi);
psi_at_centroid = nan(N_ie, N_psi);

valid_rows = cell(N_ie, 1);
triangles = cell(N_ie, 1);
edges = cell(N_ie, 1);

for idx_ie = 1:N_ie
    mask = ie_req == ie_op(idx_ie);
    rows = find(mask);
    X = [id(mask), iq(mask), psi_d(mask), ie(mask)];

    valid = all(isfinite(X), 2);
    X = X(valid, :);
    valid_rows{idx_ie} = rows(valid);

    DT = delaunayn(X(:, 1:2));
    DT_edges = get_unique_edges_from_delaunayn(DT, X(:, 1:2));
    triangles{idx_ie} = DT;
    edges{idx_ie} = DT_edges;

    for idx_psi = 1:N_psi
        points = delaunay_edge_search_N( ...
            X, 3, psi_d_values(idx_psi), DT_edges);

        contour_points{idx_ie, idx_psi} = points;
        num_points(idx_ie, idx_psi) = size(points, 1);

        if isempty(points)
            continue
        end

        mean_id(idx_ie, idx_psi) = mean(points(:, 1));
        mean_iq(idx_ie, idx_psi) = mean(points(:, 2));
        mean_ie(idx_ie, idx_psi) = mean(points(:, 4));
    end

    % A centroid does not necessarily lie on the corresponding contour.
    F = scatteredInterpolant(X(:, 1), X(:, 2), X(:, 3), ...
                             'linear', 'none');
    available = isfinite(mean_id(idx_ie, :)) & isfinite(mean_iq(idx_ie, :));
    psi_at_centroid(idx_ie, available) = F( ...
        mean_id(idx_ie, available)', mean_iq(idx_ie, available)')';
end

% Same absolute-numerator convention as Torque_LUT_Calc.m.
delta_id = diff(mean_id, 1, 1);
delta_iq = diff(mean_iq, 1, 1);
delta_ie = diff(mean_ie, 1, 1);

k_values = abs(delta_id) ./ delta_ie;
k_values(~isfinite(k_values) | abs(k_values) < eps) = NaN;
kde = mean(k_values(:, legacy_mask), 'all', 'omitnan');

k_signed = -delta_id ./ delta_ie;
k_signed(~isfinite(k_signed)) = NaN;

iq_drift = delta_iq ./ delta_ie;
iq_drift(~isfinite(iq_drift)) = NaN;

% Keep the full-range results; only kde uses the legacy fitting interval.
meta.ie_op = ie_op;
meta.ie_pairs = [ie_op(1:end-1), ie_op(2:end)];
meta.psi_d_values = psi_d_values;
meta.legacy_mask = legacy_mask;
meta.contour_points = contour_points;
meta.mean_id = mean_id;
meta.mean_iq = mean_iq;
meta.mean_ie = mean_ie;
meta.num_points = num_points;
meta.psi_at_centroid = psi_at_centroid;
meta.centroid_error = psi_at_centroid - psi_d_values;
meta.delta_id = delta_id;
meta.delta_iq = delta_iq;
meta.delta_ie = delta_ie;
meta.k_values = k_values;
meta.k_signed = k_signed;
meta.iq_drift = iq_drift;
meta.valid_rows = valid_rows;
meta.triangles = triangles;
meta.edges = edges;
end
