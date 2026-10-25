function local = investigate_local_factors(set_ids, min_rcond)
% Reuse saved 80-neighbour derivatives; no refitting or global optimization.
% Example: local = investigate_local_factors; plot_local_factor(local{1});
if nargin < 1 || isempty(set_ids)
    set_ids = 1:3;
end
if nargin < 2
    min_rcond = 1e-3;
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'matlab', 'common'));
v = load(fullfile(root, 'results', 'reduction_data.mat'), 'data');
derivatives = load(fullfile(root, 'results', 'winding_validation.mat'), 'study');
global_fit = load(fullfile(root, 'results', 'flux_map_reduction.mat'), 'study');
local = cell(size(set_ids));
thresholds = unique([1e-4 1e-3 1e-2 min_rcond]);
fprintf('80-neighbour saved derivatives; kd,kq require rcond(J) >= %g.\n', min_rcond);
fprintf('This numerical cutoff is a screening choice, not a measurement uncertainty bound.\n');

for j = 1:numel(set_ids)
    s = set_ids(j);
    [x, ~, groups, excitation] = reduction_inputs(v.data{s});
    l = derivatives.study{s}.slopes{2};
    denominator = sum(l.a.^2, 2);
    valid_local = l.valid & isfinite(l.klocal) & isfinite(denominator) & denominator > 0;
    % These are two estimates of ONE d-shift, not the two-axis pair kd,kq.
    denominator_limits = .1*median(abs(l.a), 1, 'omitnan');
    finite_ab = all(isfinite([l.a l.b]), 2);
    component_valid = finite_ab & abs(l.a) >= denominator_limits & abs(l.a) > 0 & isfinite(l.ratio);
    components = l.ratio;
    components(~component_valid) = nan;
    common = valid_local & all(component_valid, 2);
    small_Lqd = finite_ab & abs(l.a(:, 2)) < denominator_limits(2);
    condition = nan(size(groups));
    pair = nan(numel(groups), 2);
    % The two-axis solve has its own Jacobian criterion, independent of ||a||.
    finite = all(isfinite(reshape(l.L, [], 6)), 2);

    for i = find(finite)'
        J = squeeze(l.L(i, :, 1:2));
        condition(i) = rcond(J);
        if condition(i) >= min(thresholds)
            pair(i, :) = (J\l.b(i, :)')';
        end
    end
    valid_pair = finite & condition >= min_rcond & all(isfinite(pair), 2);
    f = global_fit.study{s}.p{6}(:)';
    fprintf('\nSet%d: scalar %d/%d, two-axis %d/%d; saved global F = [%.6g %.6g] A/A\n', ...
        s, nnz(valid_local), numel(groups), nnz(valid_pair), numel(groups), f);
    % The reference mask already rejects ||a|| <= 0.1 median(||a||).
    denom_limit = (.1*median(vecnorm(l.a, 2, 2), 'omitnan'))^2;
    fprintf('Saved scalar mask: a''a > %.6g (Vs/A)^2; no additional value clipping.\n', denom_limit);
    fprintf('                        min           p05        median           p95           max\n');
    names = {'k_loc', 'kd(psi_d)', 'kd(psi_q)', 'k_d (2-axis)', 'k_q (2-axis)', 'rcond(J)'};
    fields = {l.klocal(valid_local), components(component_valid(:, 1), 1), ...
        components(component_valid(:, 2), 2), pair(valid_pair, 1), pair(valid_pair, 2), condition(finite)};
    for c = 1:numel(fields)
        a = sort(fields{c});
        a = a(isfinite(a));
        q = nan(1, 5);
        if ~isempty(a)
            q = a(max(1, ceil([0 .05 .5 .95 1]*numel(a))));
            q(3) = median(a);
        end
        fprintf('%-12s %13.6g %13.6g %13.6g %13.6g %13.6g\n', names{c}, q);
    end
    fprintf('rcond cutoff       retained     median kd     median kq\n');
    for cutoff = thresholds
        keep = finite & condition >= cutoff & all(isfinite(pair), 2);
        fprintf('%12.4g %12d %13.6g %13.6g\n', cutoff, nnz(keep), median(pair(keep, :), 1, 'omitnan'));
    end
    fprintf('Median |local - saved global F|: kd %.6g, kq %.6g A/A.\n', ...
        median(abs(pair(valid_pair, :)-f), 1, 'omitnan'));
    rejected = finite & condition < min_rcond;
    if any(rejected)
        fprintf('Rejected J current ranges [id iq ie] in A: [%g %g %g] to [%g %g %g].\n', ...
            min(x(rejected, :), [], 1), max(x(rejected, :), [], 1));
    end
    pair(~valid_pair, :) = nan;
    group_summary = nan(numel(excitation), 6);
    for g = 1:numel(excitation)
        a = groups == g;
        group_summary(g, :) = [excitation(g), nnz(a & valid_local), nnz(a & valid_pair), ...
            median(l.klocal(a & valid_local), 'omitnan'), median(pair(a & valid_pair, :), 1, 'omitnan')];
    end
    disp(array2table(group_summary, 'VariableNames', ...
        {'ie_ref_A', 'scalar_points', 'pair_points', 'median_kloc', 'median_kd', 'median_kq'}));
    fprintf('Component denominator limits |Ldd| >= %.6g, |Lqd| >= %.6g Vs/A.\n', denominator_limits);
    fprintf('Component supports: psi_d %d, psi_q %d; common with LS %d/%d.\n', ...
        sum(component_valid, 1), nnz(common), numel(groups));
    fprintf('Small |Lqd|: %d points; remaining invalid psi_q factors: %d.\n', ...
        nnz(small_Lqd), nnz(~component_valid(:, 2) & ~small_Lqd));
    if any(small_Lqd)
        fprintf('Small-Lqd regions [id iq ie] A: [%g %g %g] to [%g %g %g].\n', ...
            min(x(small_Lqd, :), [], 1), max(x(small_Lqd, :), [], 1));
        fprintf('Masked |Lqd| min/median/max: %.6g %.6g %.6g Vs/A; median |iq| %.6g A.\n', ...
            min(abs(l.a(small_Lqd, 2))), median(abs(l.a(small_Lqd, 2))), ...
            max(abs(l.a(small_Lqd, 2))), median(abs(x(small_Lqd, 2))));
    end
    % All five scalar criteria below use exactly the same physical points.
    contour_k = derivatives.study{s}.k;
    factors = [components, l.klocal, repmat([contour_k l.kg], numel(groups), 1)];
    factor_names = {'psi_d'; 'psi_q'; 'local_LS'; 'contour_global'; 'derivative_global'};
    comparison = nan(5, 6);
    for c = 1:5
        a = sort(factors(common, c));
        if ~isempty(a)
            r = l.b(common, :)-factors(common, c).*l.a(common, :);
            comparison(c, :) = [numel(a), a(max(1, ceil(.05*numel(a)))), median(a), ...
                a(ceil(.95*numel(a))), median(abs(a-contour_k)), norm(r, 'fro')/norm(l.b(common, :), 'fro')];
        end
    end
    component_comparison = array2table(comparison, 'VariableNames', ...
        {'points', 'p05', 'median', 'p95', 'median_abs_delta_contour', 'relative_derivative_residual'}, ...
        'RowNames', factor_names);
    fprintf('Common-support scalar comparison (derivative residual, not flux prediction error):\n');
    disp(component_comparison);
    component_groups = nan(numel(excitation), 7);
    for g = 1:numel(excitation)
        a = groups == g;
        component_groups(g, :) = [excitation(g), nnz(a & common), nnz(a & small_Lqd), ...
            median(factors(a & common, 1:3), 1, 'omitnan'), median(l.L(a & valid_local, 1, 1), 'omitnan')];
    end
    disp(array2table(component_groups, 'VariableNames', ...
        {'ie_ref_A', 'common_points', 'small_Lqd_points', 'median_psi_d', 'median_psi_q', 'median_LS', 'median_Ldd'}));
    fprintf('Ldd/cross-coupling associations alone do not establish saturation or measurement-error causes.\n');
    local{j} = struct('set', s, 'x', x, 'groups', groups, 'excitation', excitation, ...
        'klocal', l.klocal, 'kd', pair(:, 1), 'kq', pair(:, 2), ...
        'valid_local', valid_local, 'valid_pair', valid_pair, 'rcond', condition, ...
        'min_rcond', min_rcond, 'global_F', f, 'group_summary', group_summary, ...
        'kd_from_psi_d', components(:, 1), 'kd_from_psi_q', components(:, 2), ...
        'valid_psi_d', component_valid(:, 1), 'valid_psi_q', component_valid(:, 2), ...
        'common_scalar', common, 'small_Lqd', small_Lqd, 'denominator_limits', denominator_limits, ...
        'global_contour_k', contour_k, 'global_derivative_k', l.kg, 'L', l.L, ...
        'a', l.a, 'b', l.b, 'eta', l.eta, 'residual_contour', l.residual, ...
        'residual_ls', l.b-l.klocal.*l.a, 'component_comparison', component_comparison, ...
        'component_group_summary', component_groups);
end
end
