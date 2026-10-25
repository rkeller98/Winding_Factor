function compare_flux_fits
%% Load saved models. Compare values, physical derivatives and scalar evaluation cost.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root = setup_matlab_paths('function_model');
out = fullfile(root,'results','flux_fits');
require_analysis_file(fullfile(out,'data.mat'),'prepare_flux_fit_data');
require_analysis_file(fullfile(out,'rbf_models.mat'),'fit_flux_rbf');
require_analysis_file(fullfile(out,'bspline_models.mat'),'fit_flux_bspline');
require_analysis_file(fullfile(out,'baseline_models.mat'),'fit_flux_baselines');
a = load(fullfile(out,'data.mat'));
data = a.data;
a = load(fullfile(out,'rbf_models.mat'));
rbf = a.models;
a = load(fullfile(out,'bspline_models.mat'));
spline = a.models;
a = load(fullfile(out,'baseline_models.mat'));
baseline = a.models;
labels = {'RBF','Spline','LUT','kA','kB','kC','kD','kE','kF','kG'};
study = cell(size(data));
metrics = [];
derivatives = [];
cost = [];
checks = [];
parameters = [];
lookup_checks = [];

for s = 1:numel(data)
    d = data{s};
    models = [rbf(s), spline(s), baseline(s, :)];
    if ~strcmp(d.active{3},'ie')
        models = models(1:3);
    end
    nm = numel(models);
    N = size(d.x, 1);
    pred = nan(N, 2, nm);
    J = nan(N, 2, 3, nm);
    own = false(N, nm);

    for m = 1:nm
        [pred(:, :, m), J(:, :, :, m), own(:, m)] = eval_flux_model(models{m}, d.x);
    end
    common = all(own, 2);
    masks = {common & (d.train | d.validation), common & d.test, common & d.boundary};
    eligible = {d.train | d.validation, d.test, d.boundary};

    for mode = 1:3
        for g = [0;unique(d.group)]'
            query = eligible{mode};
            if g > 0
                query = query & d.group == g;
            end
            keep = masks{mode} & query;

            for m = 1:nm
                stats = error_stats(pred(:, :, m)-d.y, keep);
                metrics = [metrics;s, mode, g, m, stats, nnz(keep)/max(1, nnz(query)), nnz(query), nnz(own(:, m) & query)];
            end
        end
    end
    % Local slopes use only final training values. They are diagnostics, not sensors.
    train = d.train | d.validation;
    ix = find(masks{2});
    q = d.x(ix, :);
    L = local_slopes(d.x(train, :), d.y(train, :), q, 80);
    L40 = local_slopes(d.x(train, :), d.y(train, :), q, 40);
    L120 = local_slopes(d.x(train, :), d.y(train, :), q, 120);
    commonL = all(isfinite(reshape(L, numel(ix), [])), 2) ...
        & all(isfinite(reshape(L40, numel(ix), [])), 2) ...
        & all(isfinite(reshape(L120, numel(ix), [])), 2);

    for m = [1:nm 11 12]
        if m <= nm
            estimate = J(ix, :, :, m);
        elseif m == 11
            estimate = L40;
        else
            estimate = L120;
        end
        for c = 1:2
            for r = 1:3
                e = 1000*(estimate(commonL, c, r)-L(commonL, c, r));
                derivatives = [derivatives;s, m, c, r, sqrt(mean(e.^2)), quantile_abs(e, .95), max_or_nan(abs(e)), nnz(commonL)];
            end
        end
    end
    % Direct centered differences check the implementation, not physical accuracy.
    checkix = ix(round(linspace(1, numel(ix), min(24, numel(ix)))));

    for m = 1:nm
        model = models{m};
        Q = d.x(checkix, :);
        [~, analytic] = eval_flux_model(model, Q);
        fd = nan(size(analytic));
        step = (model.bounds(2, :)-model.bounds(1, :))*1e-6;

        for r = 1:3
            plus = Q;
            minus = Q;
            plus(:, r) = plus(:, r)+step(r);
            minus(:, r) = minus(:, r)-step(r);
            fd(:, :, r) = (eval_flux_model(model, plus)-eval_flux_model(model, minus))/(2*step(r));
        end
        diff = fd-analytic;
        checks = [checks;s, m, max(abs(diff), [],'all','omitnan'), ...
            sqrt(mean(diff.^2,'all','omitnan')), nnz(all(isfinite(reshape(fd, size(Q, 1), [])), 2)), ...
            nnz(ismember(model.training_rows, find(d.test | d.boundary)))];
        [terms, nc, bytes] = flux_model_size(model);
        v = model;
        details = whos('v');
        parameters = [parameters;s, m, terms, nc, bytes, details.bytes];

        % All timings include two fluxes, six derivatives and the evaluator's validity work.
        queries = Q(1:min(12, size(Q, 1)), :);
        times = zeros(7, 1);

        for j = 1:30
            [yy, jj, vv] = eval_flux_model(model, queries(mod(j-1, size(queries, 1))+1, :));
        end
        for repeat = 1:7
            tic;

            for k = 1:10
                for j = 1:size(queries, 1)
                    [yy, jj, vv] = eval_flux_model(model, queries(j, :));
                end
            end
            times(repeat) = toc/(10*size(queries, 1))*1e6;
        end
        cost = [cost;s, m, median(times), max(times), min(times), 10*size(queries, 1), 7];

        if strcmp(model.kind,'lut')
            traditional = reduction_lookup(model.map, Q);
            direct = eval_flux_model(model, Q);
            lookup_checks = [lookup_checks;s, m, max(abs(direct-traditional), [],'all','omitnan')];
        end
    end
    study{s} = struct('models', {models},'prediction', pred,'J', J,'own', own,'common', common, ...
        'test_indices', ix,'slope_reference', L,'slope40', L40,'slope120', L120,'derivative_mask', commonL);
    fprintf('%s: independent test/common %d/%d; boundary supported %d/%d\n', ...
        d.name, nnz(masks{2}), nnz(d.test), nnz(masks{3}), nnz(d.boundary));
end
writetable(array2table(metrics,'VariableNames', {'dataset','mode','group','method', ...
    'rmse_d_mVs','rmse_q_mVs','p95_d_mVs','p95_q_mVs','max_d_mVs','max_q_mVs', ...
    'joint_rmse_mVs','valid_points','coverage','eligible_points','own_points'}), fullfile(out,'metrics.csv'));
writetable(array2table(derivatives,'VariableNames', {'dataset','method','output','input', ...
    'rmse_scaled','p95_scaled','max_scaled','points'}), fullfile(out,'derivatives.csv'));
writetable(array2table(cost,'VariableNames', {'dataset','method','median_us','max_batch_us', ...
    'min_batch_us','calls_per_batch','batches'}), fullfile(out,'runtime.csv'));
writetable(array2table(checks,'VariableNames', {'dataset','method','max_FD_error', ...
    'rms_FD_error','finite_queries','test_or_boundary_training_rows'}), fullfile(out,'numerical_checks.csv'));
writetable(array2table(parameters,'VariableNames', {'dataset','method','terms','flux_coefficients', ...
    'online_numeric_bytes','MATLAB_struct_bytes'}), fullfile(out,'model_sizes.csv'));
writetable(array2table(lookup_checks,'VariableNames', {'dataset','method','max_value_difference_Vs'}), ...
    fullfile(out,'lookup_consistency.csv'));
save(fullfile(out,'comparison.mat'),'study','labels','metrics','derivatives','cost','checks','parameters','-v7');
end

function L = local_slopes(x, y, q, count)
sc = max(x)-min(x);
z = x./sc;
L = nan(size(q, 1), 2, 3);

for j = 1:size(q, 1)
    [distance, ix] = mink(sum((z-q(j, :)./sc).^2, 2), min(count, size(x, 1)));
    weight = sqrt(exp(-distance/max(distance)));
    A = [ones(numel(ix), 1), z(ix, :)-q(j, :)./sc].*weight;
    [~, R] = qr(A, 0);
    if rcond(R) < 1e-7
        continue;
    end
    b = A\(y(ix, :).*weight);
    L(j, :, :) = permute(b(2:4, :)./sc', [3 2 1]);
end
end

function a = error_stats(e, keep)
r = 1000*e(keep, :);
n = size(r, 1);

if n == 0
    a = [nan(1, 7), 0];
    return;
end
a = [sqrt(mean(r.^2, 1)), quantile_abs(r(:, 1), .95), quantile_abs(r(:, 2), .95), ...
    max(abs(r), [], 1), sqrt(mean(r.^2,'all')), n];
end
function q = quantile_abs(a, p)
a = sort(abs(a(isfinite(a))));
if isempty(a)
    q = nan;
else
    q = a(ceil(p*numel(a)));
end
end
function v = max_or_nan(a)
if isempty(a)
    v = nan;
else
    v = max(a);
end
end
