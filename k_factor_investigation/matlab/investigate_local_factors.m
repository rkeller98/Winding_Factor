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
    names = {'k_loc', 'k_d', 'k_q', 'rcond(J)'};
    fields = {l.klocal(valid_local), pair(valid_pair, 1), pair(valid_pair, 2), condition(finite)};
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
    local{j} = struct('set', s, 'x', x, 'groups', groups, 'excitation', excitation, ...
        'klocal', l.klocal, 'kd', pair(:, 1), 'kq', pair(:, 2), ...
        'valid_local', valid_local, 'valid_pair', valid_pair, 'rcond', condition, ...
        'min_rcond', min_rcond, 'global_F', f, 'group_summary', group_summary);
end
end
