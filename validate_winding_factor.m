%% Local slopes, common winding factor and flux prediction (Set 1 and Set 2)
clear; clc;
root = fileparts(mfilename('fullpath'));
files = {fullfile(root, 'Resources/DataSets/Set1/extracted/evaluation/Matlab_200_201_202_FluxEval_measurements.mat'), ...
    fullfile(root, 'Resources/DataSets/Set2/extracted/evaluation/eval_flux_eesm_excitation_0A_to_40A_measurements.mat')};
neighbors = [20 40 80];
out = fullfile(root, 'results');
figures = fullfile(root, 'winding_factor_study', 'figures');
if ~isfolder(out), mkdir(out); end
if ~isfolder(figures), mkdir(figures); end
study = cell(2, 1);

for dataset = 1:2
    loaded = load(files{dataset}, 'Total'); slices = loaded.Total.pos_3000_rpm;
    names = fieldnames(slices);
    names = names(~cellfun('isempty', regexp(names, '^(pos|neg)_\d+_A$', 'once')));
    excitation = cellfun(@(name) slices.(name).Ref_exc_current, names);
    [excitation, order] = sort(excitation); names = names(order);
    x = []; psi = []; groups = [];
    for g = 1:numel(names)
        m = slices.(names{g}).Meas;
        currents = [m.PE_Id_meas_avg(:), m.PE_Iq_meas_avg(:), m.I_exc_meas(:)];
        flux = [m.PE_psi_d_meas(:), m.PE_psi_q_meas(:)];
        keep = all(isfinite([currents flux]), 2);
        x = [x; currents(keep, :)]; psi = [psi; flux(keep, :)];
        groups = [groups; g*ones(nnz(keep), 1)];
    end
    k0 = contour_ratio(x, psi, groups);
    slopes = cell(3, 1); summary = nan(3, 9);
    for h = 1:3
        L = local_slopes(x, psi, groups, neighbors(h));
        a = [L(:, 1, 1), L(:, 1, 2)]; b = [L(:, 3, 1), L(:, 3, 2)];
        an = vecnorm(a, 2, 2); bn = vecnorm(b, 2, 2);
        valid = all(isfinite([a b]), 2) & an > .1*median(an, 'omitnan');
        kLS = sum(a.*b, 2)./sum(a.^2, 2); kLS(~valid) = nan;
        kd = b(:, 1)./a(:, 1); kq = b(:, 2)./a(:, 2);
        kd(abs(a(:, 1)) < .1*median(abs(a(:, 1)), 'omitnan')) = nan;
        kq(abs(a(:, 2)) < .1*median(abs(a(:, 2)), 'omitnan')) = nan;
        eta = abs(a(:, 1).*b(:, 2)-a(:, 2).*b(:, 1))./(an.*bn);
        eta(~valid | bn < .1*median(bn, 'omitnan')) = nan;
        kGlobal = sum(a(valid, :).*b(valid, :), 'all')/sum(a(valid, :).^2, 'all');
        r0 = b-k0*a; rBest = b-kGlobal*a; rLocal = b-kLS.*a;
        sorted = sort(eta(isfinite(eta)));
        summary(h, :) = [neighbors(h), kGlobal, median(kLS, 'omitnan'), median(eta, 'omitnan'), ...
            sorted(ceil(.95*numel(sorted))), norm(r0(valid, :), 'fro')/norm(b(valid, :), 'fro'), ...
            norm(rBest(valid, :), 'fro')/norm(b(valid, :), 'fro'), ...
            norm(rLocal(valid, :), 'fro')/norm(b(valid, :), 'fro'), nnz(valid)];
        slopes{h} = struct('L', L, 'kd', kd, 'kq', kq, 'kLS', kLS, 'eta', eta, 'r0', r0, 'valid', valid);
        fprintf('Set %g, %g neighbors: k0 %.6f, global k %.6f, relative residual %.4f\n', ...
            dataset, neighbors(h), k0, kGlobal, summary(h, 6));
    end

    % Compare neighborhood sizes on the same current points.
    derivative_common = slopes{1}.valid & slopes{2}.valid & slopes{3}.valid;
    sensitivity = nan(3, 3);
    for h = 1:3
        L = slopes{h}.L;
        a = [L(derivative_common, 1, 1), L(derivative_common, 1, 2)];
        b = [L(derivative_common, 3, 1), L(derivative_common, 3, 2)];
        k = sum(a.*b, 'all')/sum(a.^2, 'all');
        sensitivity(h, :) = [k, median(slopes{h}.kLS(derivative_common)), norm(b-k*a, 'fro')/norm(b, 'fro')];
    end

    %% Local maps at the nominal 20 A slice; no additional gridding
    fig = figure(10+dataset); clf; theme(fig, 'light');
    set(fig, 'Position', [50 50 1300 760]); t = tiledlayout(2, 3, 'TileSpacing', 'compact');
    g20 = find(excitation == 20); s = groups == g20; local = slopes{2};
    values = {local.kd, local.kq, local.kLS, local.eta, 1000*vecnorm(local.r0, 2, 2), 1000*local.L(:, 1, 2)};
    titles = {'k_d^{(\psi_d)} in A/A', 'k_d^{(\psi_q)} in A/A', 'Local k_{LS} in A/A', ...
        'Normalized determinant in 1', 'Constant-k residual norm in mH', 'L_{qd} in mH'};
    for j = 1:6
        nexttile; scatter(x(s, 1), x(s, 2), 15, values{j}(s), 'filled');
        xlabel('I_d in A'); ylabel('I_q in A'); title(titles{j}); colorbar; grid on;
    end
    title(t, sprintf('Set %g: nominal excitation current 20 in A, 40 neighbors', dataset));
    exportgraphics(t, fullfile(figures, sprintf('local_set%g.png', dataset)), 'Resolution', 180);

    %% Hold out one entire interior excitation group at a time
    prediction = nan(size(psi, 1), 2, 3); fold_k = nan(numel(names), 2);
    for g = 2:numel(names)-1
        train = groups ~= g; test = groups == g;
        xt = x(train, :); yt = psi(train, :); gt = groups(train);
        L = local_slopes(xt, yt, gt, 40);
        a = [L(:, 1, 1), L(:, 1, 2)]; b = [L(:, 3, 1), L(:, 3, 2)];
        an = vecnorm(a, 2, 2); valid = all(isfinite([a b]), 2) & an > .1*median(an, 'omitnan');
        fold_k(g, :) = [contour_ratio(xt, yt, gt), sum(a(valid, :).*b(valid, :), 'all')/sum(a(valid, :).^2, 'all')];
        for model = 1:3
            if model < 3
                k = fold_k(g, model);
                points = [xt(:, 1)+k*xt(:, 3), xt(:, 2)];
                query = [x(test, 1)+k*x(test, 3), x(test, 2)];
            else
                scale = max(xt)-min(xt); points = xt./scale; query = x(test, :)./scale;
            end
            F = scatteredInterpolant(points, yt, 'linear', 'none');
            prediction(test, :, model) = F(query);
        end
        fprintf('Set %g: held-out excitation %g in A\n', dataset, excitation(g));
    end
    interior = groups > 1 & groups < numel(names);
    common = interior & all(isfinite(prediction), [2 3]);
    cv = nan(3, 5); fold_rmse = nan(numel(names), 2, 3);
    for model = 1:3
        error = 1000*(prediction(:, :, model)-psi);
        cv(model, :) = [sqrt(mean(error(common, :).^2)), max(abs(error(common, :))), nnz(common)];
        for g = 2:numel(names)-1
            use = common & groups == g;
            fold_rmse(g, :, model) = sqrt(mean(error(use, :).^2));
        end
    end
    fig = figure(20+dataset); clf; theme(fig, 'light');
    set(fig, 'Position', [50 50 1100 450]); t = tiledlayout(1, 2, 'TileSpacing', 'compact');
    labels = {'Reduced, contour k', 'Reduced, global LS k', 'Full 3D interpolation'};
    for component = 1:2
        nexttile; plot(excitation, squeeze(fold_rmse(:, component, :)), '.-', 'LineWidth', 1);
        xlabel('Held-out excitation current in A'); ylabel('Flux RMSE in mWb');
        title(sprintf('Set %g, %s flux', dataset, char('d'+13*(component-1)))); legend(labels, 'Location', 'best'); grid on;
    end
    exportgraphics(t, fullfile(figures, sprintf('prediction_set%g.png', dataset)), 'Resolution', 180);

    %% Grid roundtrip: same all-data scattered field, different raster sizes
    projected = [x(:, 1)+k0*x(:, 3), x(:, 2)];
    F = scatteredInterpolant(projected, psi, 'linear', 'none');
    grid_sizes = [101 201 401]; returned = nan(size(psi, 1), 2, 3);
    for h = 1:3
        d = linspace(min(projected(:, 1)), max(projected(:, 1)), grid_sizes(h));
        q = linspace(min(projected(:, 2)), max(projected(:, 2)), grid_sizes(h));
        [D, Q] = meshgrid(d, q); grid_psi = F([D(:), Q(:)]);
        for component = 1:2
            values = reshape(grid_psi(:, component), size(D));
            returned(:, component, h) = interp2(d, q, values, projected(:, 1), projected(:, 2), 'linear', nan);
        end
    end
    supported = all(isfinite(returned), [2 3]); raster = nan(3, 5);
    for h = 1:3
        error = 1000*(returned(supported, :, h)-psi(supported, :));
        raster(h, :) = [sqrt(mean(error.^2)), max(abs(error)), nnz(supported)];
    end
    study{dataset} = struct('x', x, 'psi', psi, 'groups', groups, 'excitation', excitation, ...
        'k0', k0, 'slopes', {slopes}, 'summary', summary, 'prediction', prediction, ...
        'fold_k', fold_k, 'common', common, 'cv', cv, 'fold_rmse', fold_rmse, 'raster', raster, ...
        'sensitivity', sensitivity, 'derivative_common', derivative_common);
end

%% Compact comparisons for the manuscript
fig = figure(30); clf; theme(fig, 'light'); set(fig, 'Position', [50 50 1100 450]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact');
for dataset = 1:2
    nexttile; plot(neighbors, study{dataset}.sensitivity(:, 1:2), '.-', 'LineWidth', 1); hold on;
    yline(study{dataset}.k0, '--'); xlabel('Number of neighbors'); ylabel('k in A/A');
    title(sprintf('Set %g', dataset)); legend('Global LS', 'Median local LS', 'MeasEval', 'Location', 'best'); grid on;
end
exportgraphics(t, fullfile(figures, 'sensitivity.png'), 'Resolution', 180);
fig = figure(31); clf; theme(fig, 'light'); set(fig, 'Position', [50 50 1100 450]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact');
for dataset = 1:2
    nexttile; plot(grid_sizes, study{dataset}.raster(:, 1:2), '.-', 'LineWidth', 1);
    xlabel('Grid points per axis'); ylabel('Roundtrip RMSE in mWb'); title(sprintf('Set %g', dataset));
    legend('d flux', 'q flux'); grid on;
end
exportgraphics(t, fullfile(figures, 'raster.png'), 'Resolution', 180);
save(fullfile(out, 'winding_validation.mat'), 'study', 'neighbors', 'files');
disp('Column order: neighbors, global k, median local k, median eta, p95 eta, residual k0/global/local, valid points');
disp(study{1}.summary); disp(study{2}.summary);
disp('Prediction: d/q RMSE in mWb, d/q maximum in mWb, common supported points');
disp(study{1}.cv); disp(study{2}.cv);

function L = local_slopes(x, psi, groups, count)
% Fit both flux components independently: psi = c + Ld*dId + Lq*dIq + Le*dIe.
    scale = max(x)-min(x); z = x./scale; L = nan(size(x, 1), 3, 2);
    for p = 1:size(x, 1)
        [distance, index] = mink(sum((z-z(p, :)).^2, 2), count);
        if numel(unique(groups(index))) < 2, continue; end
        weight = sqrt(exp(-distance/max(distance)));
        A = [ones(count, 1), z(index, :)-z(p, :)].*weight;
        [~, R] = qr(A, 0);
        if rcond(R) < 1e-6, continue; end
        coefficients = A\(psi(index, :).*weight);
        L(p, :, :) = reshape(coefficients(2:4, :)./scale', [1 3 2]);
    end
end

function k = contour_ratio(x, psi, groups)
% Original MeasEval contour-centroid calculation, including its absolute value.
    refs = unique(groups); means = nan(numel(refs), 100, 3);
    low = min(psi(:, 1)); span = max(psi(:, 1))-low;
    levels = linspace(low+.25*span, low+.75*span, 100);
    for g = 1:numel(refs)
        use = groups == refs(g); c = x(use, :); f = psi(use, 1);
        tri = delaunay(c(:, 1), c(:, 2));
        edge = unique(sort([tri(:, [1 2]); tri(:, [2 3]); tri(:, [1 3])], 2), 'rows');
        length2 = sum((c(edge(:, 1), 1:2)-c(edge(:, 2), 1:2)).^2, 2);
        edge = edge(length2 < 9*mean(length2), :); a = edge(:, 1); b = edge(:, 2);
        for p = 1:100
            t = (levels(p)-f(a))./(f(b)-f(a)); hit = t >= 0 & t <= 1;
            means(g, p, :) = mean(c(a(hit), :)+t(hit).*(c(b(hit), :)-c(a(hit), :)), 1);
        end
    end
    ratio = abs(diff(means(:, :, 1)))./diff(means(:, :, 3));
    ratio(~isfinite(ratio) | ratio < eps) = nan; k = mean(ratio(:), 'omitnan');
end
