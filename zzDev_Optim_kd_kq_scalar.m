eval_set = load_data();
meas = eval_set.measurement;

%% Analytical pointwise inverse of the reference excitation slice
id = meas.id_req(:);
iq = meas.iq_req(:);
ie = meas.ie_req(:);
psid = meas.psi_d(:);
psiq = meas.psi_q(:);

valid = isfinite(id) & isfinite(iq) & isfinite(ie) & ...
        isfinite(psid) & isfinite(psiq);

id = id(valid);
iq = iq(valid);
ie = ie(valid);
psid = psid(valid);
psiq = psiq(valid);

if isempty(ie)
    error('No finite measurement points are available.');
end

ie_vec = unique(ie);
ie_ref = min(ie_vec);
delta_ie = ie - ie_ref;
reference = ie == ie_ref;
P = [id(reference), iq(reference)];
F = [psid(reference), psiq(reference)];

% Repeated current coordinates can define a surface only if their fluxes
% agree to roundoff. Keep the observations; never average conflicting data.
[P_ref, first, group] = unique(P, 'rows');
counts = accumarray(group, 1);
F_min = [accumarray(group, F(:, 1), [], @min), ...
         accumarray(group, F(:, 2), [], @min)];
F_max = [accumarray(group, F(:, 1), [], @max), ...
         accumarray(group, F(:, 2), [], @max)];
flux_roundoff = 32 * eps(max(abs(F), [], 1));
conflicting_group = any(F_max - F_min > flux_roundoff, 2);
duplicate_reference = table(P_ref(:, 1), P_ref(:, 2), counts, ...
    F_min(:, 1), F_max(:, 1), F_min(:, 2), F_max(:, 2), ...
    conflicting_group, 'VariableNames', {'id', 'iq', 'Count', ...
    'PsiDMin', 'PsiDMax', 'PsiQMin', 'PsiQMax', 'Conflicting'});
duplicate_reference = duplicate_reference(counts > 1, :);
fprintf('Reference excitation: %g A; %d points, %d unique coordinates.\n', ...
    ie_ref, size(P, 1), size(P_ref, 1));
if ~isempty(duplicate_reference)
    disp(duplicate_reference);
end
if any(conflicting_group)
    error(['%d reference coordinate groups contain conflicting fluxes. ' ...
        'Choose the measurement condition explicitly before identification; ' ...
        'see duplicate_reference. No observations have been averaged.'], ...
        nnz(conflicting_group));
end
F_ref = F(first, :);
fprintf('Collapsed %d roundoff-equivalent reference repeats for triangulation.\n', ...
    size(P, 1) - size(P_ref, 1));

if size(P_ref, 1) < 3
    error('At least three distinct reference current coordinates are required.');
end
reference_mesh = delaunayTriangulation(P_ref);
tri = reference_mesh.ConnectivityList;
if isempty(tri)
    error('The reference current coordinates do not span a two-dimensional domain.');
end

[S, psi_rec, triangle_id, barycentric, ambiguous, degenerate, outside, ...
    excluded_triangles] = invert_reference(P_ref, F_ref, tri, [psid, psiq]);

sd = S(:, 1);
sq = S(:, 2);
psid_rec = psi_rec(:, 1);
psiq_rec = psi_rec(:, 2);
valid_inverse = isfinite(sd) & isfinite(sq);
kd = nan(size(id));
kq = nan(size(iq));
valid_factor = valid_inverse & ~reference;
kd(valid_factor) = (sd(valid_factor) - id(valid_factor)) ./ delta_ie(valid_factor);
kq(valid_factor) = (sq(valid_factor) - iq(valid_factor)) ./ delta_ie(valid_factor);
factor_nonfinite = valid_factor & (~isfinite(kd) | ~isfinite(kq));
kd(factor_nonfinite) = NaN;
kq(factor_nonfinite) = NaN;
valid_factor(factor_nonfinite) = false;
residual_d = psid_rec - psid;
residual_q = psiq_rec - psiq;

% All rows retain their original position in the finite-filtered measurements.
results = table(id, iq, ie, kd, kq, sd, sq, psid, psiq, psid_rec, psiq_rec, ...
    residual_d, residual_q, reference, valid_inverse, valid_factor, ...
    outside, degenerate, ambiguous, factor_nonfinite, triangle_id, barycentric);

%% Slice diagnostics (errors describe accepted inverses only)
n_slices = numel(ie_vec);
input_count = zeros(n_slices, 1);
inverse_count = input_count;
factor_count = input_count;
outside_count = input_count;
degenerate_count = input_count;
ambiguous_count = input_count;
nonfinite_count = input_count;
rmse_d = nan(n_slices, 1);
rmse_q = rmse_d;
max_abs_d = rmse_d;
max_abs_q = rmse_d;
max_abs_kd = rmse_d;
max_abs_kq = rmse_d;
for j = 1:n_slices
    slice = ie == ie_vec(j);
    accepted = slice & valid_inverse;
    factors = slice & valid_factor;
    input_count(j) = nnz(slice);
    inverse_count(j) = nnz(accepted);
    factor_count(j) = nnz(factors);
    outside_count(j) = nnz(slice & outside);
    degenerate_count(j) = nnz(slice & degenerate);
    ambiguous_count(j) = nnz(slice & ambiguous);
    nonfinite_count(j) = nnz(slice & factor_nonfinite);
    if any(accepted)
        rmse_d(j) = sqrt(mean(residual_d(accepted).^2));
        rmse_q(j) = sqrt(mean(residual_q(accepted).^2));
        max_abs_d(j) = max(abs(residual_d(accepted)));
        max_abs_q(j) = max(abs(residual_q(accepted)));
    end
    if any(factors)
        max_abs_kd(j) = max(abs(kd(factors)));
        max_abs_kq(j) = max(abs(kq(factors)));
    end
end
slice_summary = table(ie_vec, input_count, inverse_count, factor_count, ...
    outside_count, degenerate_count, ambiguous_count, nonfinite_count, ...
    rmse_d, max_abs_d, rmse_q, max_abs_q, max_abs_kd, max_abs_kq);
fprintf('Excluded %d of %d flux triangles as singular or ill-conditioned.\n', ...
    nnz(excluded_triangles), size(tri, 1));
fprintf('Residual RMSE and maximum absolute errors are in Vs.\n');
fprintf('Degenerate and ambiguous masks can overlap; neither is accepted.\n');
disp(slice_summary);

%% Scatter plots use physical current coordinates; no factor interpolation
plot_slices(id, iq, ie, ie_vec, [kd, kq], ...
    {'k_d', 'k_q'}, {'k_d in A/A', 'k_q in A/A'});
plot_slices(id, iq, ie, ie_vec, [residual_d, residual_q], ...
    {'d flux reconstruction residual', 'q flux reconstruction residual'}, ...
    {'Flux linkage residual in Vs', 'Flux linkage residual in Vs'});


function [S, reconstructed, triangle_id, barycentric, ambiguous, ...
        degenerate, outside, excluded] = invert_reference(P, F, tri, target)
    n = size(target, 1);
    S = nan(n, 2);
    reconstructed = nan(n, 2);
    triangle_id = nan(n, 1);
    barycentric = nan(n, 3);
    ambiguous = false(n, 1);
    degenerate = false(n, 1);
    excluded = false(size(tri, 1), 1);
    candidate_seen = false(n, 1);
    coordinate_error = zeros(n, 2);

    % Normalize each flux component before conditioning/membership checks.
    % 32 eps allows accumulated roundoff, not a physical measurement tolerance.
    % Require half the digits both in conditioning and propagated roundoff.
    flux_scale = max(abs(F), [], 1);
    flux_scale(flux_scale == 0) = 1;
    V = F ./ flux_scale;
    Y = target ./ flux_scale;
    roundoff = 32 * eps;
    current_roundoff = 32 * eps(max(abs(P), [], 1));

    for t = 1:size(tri, 1)
        vertices = tri(t, :);
        Ft = F(vertices, :);
        Pt = P(vertices, :);
        Vt = V(vertices, :);
        % Form differences before scaling to avoid subtracting rounded ratios.
        B = ((Ft(2:3, :) - Ft(1, :)) ./ flux_scale).';
        quality = rcond(B);
        excluded(t) = ~isfinite(quality) || quality <= sqrt(eps);
        weight_tolerance = Inf;
        if ~excluded(t)
            % B\eye(2) is a linear solve. Account for small triangles relative
            % to the flux scale as well as their shape; rcond alone misses this.
            weight_tolerance = roundoff * (1 + norm(B \ eye(2), inf));
            excluded(t) = weight_tolerance >= sqrt(eps);
        end
        nearby = find(all(Y >= min(Vt, [], 1) - roundoff & ...
                          Y <= max(Vt, [], 1) + roundoff, 2));
        if isempty(nearby)
            continue
        end
        if excluded(t)
            % Preserve the support of rejected triangles in diagnostics, even
            % for a collapsed triangle (a segment or a point in flux space).
            Z = Y(nearby, :);
            touching = inpolygon(Z(:, 1), Z(:, 2), Vt(:, 1), Vt(:, 2));
            for edge = 1:3
                a = Vt(edge, :);
                b = Vt(mod(edge, 3) + 1, :);
                v = b - a;
                if dot(v, v) == 0
                    distance = vecnorm(Z - a, 2, 2);
                else
                    fraction = max(0, min(1, (Z - a) * v.' / dot(v, v)));
                    distance = vecnorm(Z - (a + fraction * v), 2, 2);
                end
                touching = touching | distance <= roundoff;
            end
            degenerate(nearby(touching)) = true;
            continue
        end

        rhs = ((target(nearby, :) - Ft(1, :)) ./ flux_scale).';
        lambda23 = (B \ rhs).';
        weights = [1 - sum(lambda23, 2), lambda23];
        inside = all(weights >= -weight_tolerance & ...
                     weights <= 1 + weight_tolerance, 2);
        idx = nearby(inside);
        weights = weights(inside, :);
        if isempty(idx)
            continue
        end
        candidate = weights * Pt;
        error_bound = weight_tolerance * sum(abs(Pt(2:3, :) - Pt(1, :)), 1) ...
            + current_roundoff;
        seen = candidate_seen(idx);
        % Adjacent triangles can give the same inverse on shared edges. Only
        % differences exceeding the propagated current roundoff are ambiguous.
        different = seen & any(abs(candidate - S(idx, :)) > ...
            coordinate_error(idx, :) + error_bound, 2);
        ambiguous(idx(different)) = true;
        fresh = idx(~seen);
        S(fresh, :) = candidate(~seen, :);
        reconstructed(fresh, :) = weights(~seen, :) * Ft;
        triangle_id(fresh) = t;
        barycentric(fresh, :) = weights(~seen, :);
        coordinate_error(fresh, :) = repmat(error_bound, numel(fresh), 1);
        candidate_seen(idx) = true;
    end

    outside = ~candidate_seen & ~degenerate;
    unavailable = ambiguous | degenerate | outside;
    S(unavailable, :) = NaN;
    reconstructed(unavailable, :) = NaN;
    triangle_id(unavailable) = NaN;
    barycentric(unavailable, :) = NaN;
end


function plot_slices(id, iq, ie, ie_vec, values, names, labels)
    for c = 1:size(values, 2)
        fig = figure('Name', names{c});
        tabs = uitabgroup(fig);
        finite_values = values(isfinite(values(:, c)), c);
        limits = [min(finite_values), max(finite_values)];
        for j = 1:numel(ie_vec)
            tab = uitab(tabs, 'Title', sprintf('i_e = %g A', ie_vec(j)));
            ax = axes(tab);
            slice = ie == ie_vec(j);
            available = slice & isfinite(values(:, c));
            scatter(ax, id(slice), iq(slice), 12, [0.75, 0.75, 0.75], 'x');
            hold(ax, 'on');
            scatter(ax, id(available), iq(available), 28, values(available, c), 'filled');
            xlabel(ax, 'd current in A');
            ylabel(ax, 'q current in A');
            title(ax, sprintf('%s, i_e = %g A (%d/%d available; gray = unavailable)', ...
                names{c}, ie_vec(j), nnz(available), nnz(slice)));
            grid(ax, 'on');
            axis(ax, 'equal');
            colormap(ax, turbo(256));
            cb = colorbar(ax);
            cb.Label.String = labels{c};
            if numel(limits) == 2 && limits(1) < limits(2)
                clim(ax, limits);
            end
        end
    end
end
