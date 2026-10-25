function study = analyze_flux_map_reduction(dataset_selection)
% A--H: controlled reconstruction and nested whole-excitation-group prediction.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root = setup_matlab_paths('k_factor');
out = fullfile(root,'results');
require_analysis_file(fullfile(out,'reduction_data.mat'),'run_matlab_analysis([],[],''prepare'')');
loaded = load(fullfile(out,'reduction_data.mat'),'data');
data = loaded.data;

if nargin == 0
    dataset_selection = 1:numel(data);
end
labels = {'A contour','B fixed-q d fit','C joint constant','D excitation shift', ...
    'E q/excitation shift','F two-axis constant','G three-current shift','H full three-current'};
n_grid = 81;
study = cell(size(data));
parameters = [1 1 1 2 3 2 5 0];

for s = dataset_selection
    [x, y, groups, excitation] = reduction_inputs(data{s});
    ng = numel(excitation);
    fprintf('\nSET %d: all-data parameter fitting\n', s);
    tic;
    [p, fitinfo, current_scale] = fit_methods(x, y, groups, n_grid);
    maps = cell(8, 1);
    reconstruction = nan(size(y, 1), 2, 8);

    for m = 1:8
        maps{m} = fit_reduction_map(x, y, m, p{m}, current_scale, n_grid, true);
        reconstruction(:, :, m) = reduction_lookup(maps{m}, x);
    end
    prediction = nan(size(reconstruction));
    fold_p = cell(ng, 8);
    fold_info = cell(ng, 1);

    for g = 1:ng
        train = groups ~= g;
        test = ~train;
        fprintf('Set%d: outer holdout %g A (%d/%d)\n', s, excitation(g), g, ng);
        [pg, info, sc] = fit_methods(x(train, :), y(train, :), groups(train), n_grid);

        for m = 1:8
            map = fit_reduction_map(x(train, :), y(train, :), m, pg{m}, sc, n_grid, true);
            prediction(test, :, m) = reduction_lookup(map, x(test, :));
            fold_p{g, m} = pg{m};
        end
        fold_info{g} = info;
        save(fullfile(out, sprintf('Set%d_progress.mat', s)),'p','fold_p','prediction','g','fitinfo','fold_info');
    end
    recon_common = all(isfinite(reconstruction), [2 3]);
    inner = groups > 1 & groups < ng;
    cv_common = all(isfinite(prediction), [2 3]) & inner;
    boundary_common = all(isfinite(prediction(:, :, 1:7)), [2 3]) & ~inner;
    metrics = [];
    pergroup = [];

    for mode = 1:2
        pred = reconstruction;
        common = recon_common;
        eligible = true(size(groups));

        if mode == 2
            pred = prediction;
            common = cv_common;
            eligible = inner;
        end
        for m = 1:8
            metrics = [metrics;s, mode, m, stats(pred(:, :, m)-y, common, nnz(eligible)), ...
                nnz(eligible & all(isfinite(pred(:, :, m)), 2)), parameters(m), ...
                maps{m}.nodes, 16*maps{m}.nodes+8*numel(p{m})];

            for g = 1:ng
                pergroup = [pergroup;s, mode, m, excitation(g), ...
                    stats(pred(:, :, m)-y, common & groups == g, nnz(groups == g)), ...
                    nnz(groups == g & all(isfinite(pred(:, :, m)), 2))];
            end
        end
    end
    boundary = [];

    for m = 1:8
        common = boundary_common & all(isfinite(prediction(:, :, m)), 2);
        boundary = [boundary;s, m, stats(prediction(:, :, m)-y, common, nnz(~inner)), ...
            stats(prediction(:, :, m)-y, ~inner & all(isfinite(prediction(:, :, m)), 2), nnz(~inner))];
    end
    scaling = nan(size(y, 1), 2, 4);

    for j = 1:4
        m = 1;
        if j > 2
            m = 8;
        end
        map = fit_reduction_map(x, y, m, p{m}, current_scale, n_grid, mod(j, 2) == 1);
        scaling(:, :, j) = reduction_lookup(map, x);
    end
    scale_common = all(isfinite(scaling), [2 3]);
    scaling_stats = [];

    for j = 1:4
        scaling_stats = [scaling_stats;j, stats(scaling(:, :, j)-y, scale_common, size(x, 1))];
    end
    raster_stats = [];
    raster = nan(size(y, 1), 2, 3);
    sizes = [41 81 161];

    for j = 1:3
        map = fit_reduction_map(x, y, 1, p{1}, current_scale, sizes(j), true);
        raster(:, :, j) = reduction_lookup(map, x);
    end
    raster_common = all(isfinite(raster), [2 3]);

    for j = 1:3
        raster_stats = [raster_stats;sizes(j), stats(raster(:, :, j)-y, raster_common, size(x, 1))];
    end
    jacobian = 1+x(:, 3).*(p{7}(4)/current_scale(1)+2*p{7}(5)*x(:, 1)/current_scale(1)^2);
    xref = x;
    xref(:, 3) = excitation(groups);
    reference_prediction = reduction_lookup(maps{1}, xref);
    refmask = recon_common & all(isfinite(reference_prediction), 2);
    reference_stats = stats(reference_prediction-reconstruction(:, :, 1), refmask, size(x, 1));
    result = struct('p', {p},'fitinfo', {fitinfo},'maps', {maps},'current_scale', current_scale, ...
        'reconstruction', reconstruction,'prediction', prediction,'fold_p', {fold_p}, ...
        'fold_info', {fold_info},'recon_common', recon_common,'cv_common', cv_common, ...
        'boundary_common', boundary_common,'metrics', metrics,'pergroup', pergroup, ...
        'boundary', boundary,'scaling_stats', scaling_stats,'raster_stats', raster_stats, ...
        'jacobian', jacobian,'reference_stats', reference_stats,'elapsed_seconds', toc,'fit_revision', 2);
    study{s} = result;
    save(fullfile(out, sprintf('Set%d_reduction.mat', s)),'result','labels','parameters','n_grid','-v7');
    fprintf('Set%d completed in %.1f seconds\n', s, result.elapsed_seconds);
end
% Individual sets may be executed concurrently; collect only finished new runs.
for s = 1:numel(data)
    file = fullfile(out, sprintf('Set%d_reduction.mat', s));

    if ~isfile(file)
        return;
    end
    v = load(file,'result');

    if ~isscalar(v.result) || ~isfield(v.result,'fit_revision')
        return;
    end
    study{s} = v.result;
end
save(fullfile(out,'flux_map_reduction.mat'),'study','labels','parameters','n_grid','-v7');
allmetrics = cellfun(@(r) r.metrics, study,'UniformOutput', false);
allmetrics = vertcat(allmetrics{:});
allgroups = cellfun(@(r) r.pergroup, study,'UniformOutput', false);
allgroups = vertcat(allgroups{:});
names = {'set','mode','method','rmse_d_mVs','rmse_q_mVs','joint_rmse_mVs', ...
    'max_abs_mVs','p95_vector_mVs','valid_points','coverage','own_points','shift_parameters','grid_nodes','flux_bytes'};
writetable(array2table(allmetrics,'VariableNames', names), fullfile(out,'method_metrics.csv'));
writetable(array2table(allgroups,'VariableNames', [names(1:3), {'ie_nominal_A'}, names(4:11)]), fullfile(out,'group_metrics.csv'));
end

function [p, info, sc] = fit_methods(x, y, g, n)
sc = max(abs(x), [], 1);
k = estimate_contour_factor(x, y, g);
unit = max(abs(k), 1);
ug = unique(g);
inner = ug(2:end-1);
folds = unique(inner(round(linspace(1, numel(inner), min(3, numel(inner))))));
parts = cell(numel(folds), 1);

for j = 1:numel(folds)
    tr = g ~= folds(j);
    ix = find(g == folds(j));
    ix = ix(round(linspace(1, numel(ix), min(120, numel(ix)))));
    base = fit_reduction_map(x(tr, :), y(tr, :), 1, k, sc, n, true);
    keep = all(isfinite(reduction_lookup(base, x(ix, :))), 2);
    ix = ix(keep);
    parts{j} = struct('x', x(tr, :),'y', y(tr, :),'query', x(ix, :),'target', y(ix, :));
end
p = cell(8, 1);
p{1} = k;
p{8} = [];
info = cell(8, 1);
info{1} = struct('inner_groups', folds);

for m = 2:3
    fun = @(v) objective(v*unit, m, parts, sc, n, k);
    seeds = linspace(-.5, 2.5, 9);
    score = arrayfun(fun, seeds);
    [~, order] = sort(score);
    best = inf;
    chosen = nan;
    total = 0;

    for j = 1:2
        h = order(j);
        lo = seeds(max(1, h-1));
        hi = seeds(min(numel(seeds), h+1));
        [v, fval, exitflag, output] = fminbnd(fun, lo, hi, optimset('TolX', .002,'Display','off'));
        total = total+output.funcCount;

        if fval < best
            best = fval;
            chosen = v;
            flag = exitflag;
        end
    end
    p{m} = chosen*unit;
    info{m} = struct('objective_Vs2', best,'exitflag', flag, ...
        'evaluations', total+9,'inner_groups', folds,'search_range', [-.5 2.5]*unit);
end
count = [1 1 1 2 3 2 5];

for m = 4:7
    base = zeros(1, count(m));
    base(1) = p{3}/unit;

    % Nonzero search coordinates give every coefficient the same initial step.
    % A zero-valued physical coefficient must not get MATLAB's tiny zero step.
    start = ones(1, count(m));
    fun = @(v) objective((base+v-1)*unit, m, parts, sc, n, k);
    [v, fval, flag, output] = fminsearch(fun, start, ...
        optimset('MaxIter', 220,'MaxFunEvals', 350,'TolX', .003,'TolFun', 1e-10,'Display','off'));
    p{m} = (base+v-1)*unit;
    info{m} = struct('objective_Vs2', fval,'exitflag', flag, ...
        'evaluations', output.funcCount,'inner_groups', folds);
end
end

function cost = objective(p, m, parts, sc, n, k)
unit = max(abs(k), 1);
v = p/unit;

if any(abs(v) > 4)
    cost = 1+sum(v.^2);
    return;
end
if m == 7
    [D, E] = ndgrid([-sc(1) sc(1)], [-sc(3) sc(3)]);
    jac = 1+E(:).*(p(4)/sc(1)+2*p(5)*D(:)/sc(1)^2);

    if min(jac) < .25
        cost = 1+sum((jac-.25).^2);
        return;
    end
end
cost = 0;

for j = 1:numel(parts)
    d = parts{j};
    map = fit_reduction_map(d.x, d.y, m, p, sc, n, true);
    pred = reduction_lookup(map, d.query);
    r = pred-d.target;

    if m == 2
        r = r(:, 1);
    end
    r(~isfinite(r)) = 1;
    % Fixed denominator: lost support cannot improve the fit.
    cost = cost+mean(r(:).^2)/numel(parts);
end
if m >= 4
    cost = cost+1e-8*sum(v(2:end).^2);
end
end

function a = stats(r, mask, total)
r = 1000*r(mask, :);
n = size(r, 1);

if n == 0
    a = [nan nan nan nan nan 0 0];
    return;
end
v = sort(vecnorm(r, 2, 2));
a = [sqrt(mean(r.^2, 1)), sqrt(mean(r.^2,'all')), ...
    max(abs(r), [],'all'), v(ceil(.95*n)), n, n/max(total, 1)];
end
