%% Saved MeasEval slices versus the merged two-dimensional flux map
% Run from any directory. Uses actual measured currents; no extrapolation.
root = fileparts(mfilename('fullpath'));
files = {fullfile(root, 'Resources/DataSets/Set1/FluxEval/Matlab_200_201_202_FluxEval.mat'), ...
    fullfile(root, 'Resources/DataSets/Set2/FluxEval/eval_flux_eesm_excitation_0A_to_40A.mat')};
out = fullfile(root, 'results'); if ~isfolder(out), mkdir(out); end
figures = fullfile(root, 'winding_factor_study/figures');
if ~isfolder(figures), mkdir(figures); end
diagnosis = cell(1, 2);

for dataset = 1:2
    data = load(files{dataset}, 'Total', 'LUT');
    map = data.LUT.CtrlMaps; complete = map.complete;
    k0 = map.EESM_winding_ratio; blocks = data.Total.pos_3000_rpm;
    names = fieldnames(blocks);
    refs = cellfun(@(name) blocks.(name).Ref_exc_current, names);
    [excitation, order] = sort(refs); names = names(order);
    x = []; psi = []; returned = []; formula = []; groups = [];
    coordinate = strings(numel(names), 1);
    axis_error = nan(numel(names), 2);

    %% A: measured flux. B: lookup of its own excitation-reference slice.
    for g = 1:numel(names)
        block = blocks.(names{g}); m = block.Meas;
        current = [m.PE_Id_meas_avg(:), m.PE_Iq_meas_avg(:), m.I_exc_meas(:)];
        flux = [m.PE_psi_d_meas(:), m.PE_psi_q_meas(:)];
        keep = all(isfinite([current flux]), 2);
        slice = find(abs(complete.I_exc_vec-excitation(g)) < 1e-6);
        d = complete.Id_breakpoints(slice, :); q = complete.Iq_breakpoints(slice, :);
        physical = linspace(min(block.Calc.Grid.Id_fine), max(block.Calc.Grid.Id_fine), numel(d));
        offset = complete.EESM_winding_ratio*excitation(g);
        axis_error(g, :) = [max(abs(d-physical)), max(abs(d-physical-offset))];
        % New metadata is authoritative. Older files are identified by axis alignment.
        if isfield(complete, 'Coordinates')
            shifted = strcmp(complete.Coordinates.d_axis, 'Id_star');
        else
            [distance, convention] = min(axis_error(g, :));
            if distance > 1e-6, error('Unknown complete-map coordinates in %s', names{g}); end
            shifted = convention == 2;
        end
        query_d = current(:, 1);
        coordinate(g) = "physical Id";
        if shifted
            % This is a relabelled physical slice: use its fixed REFERENCE offset.
            query_d = query_d+offset; coordinate(g) = "Id_star, slice reference";
        end
        values = nan(size(flux));
        values(:, 1) = interp2(d, q, complete.Psi_d(:, :, slice), query_d, current(:, 2), 'linear', nan);
        values(:, 2) = interp2(d, q, complete.Psi_q(:, :, slice), query_d, current(:, 2), 'linear', nan);
        source_formula = [(m.PE_Uq_meas_avg(:)-m.Const_Machine_R_s(:).*current(:, 2))./m.PE_w_el_rad_avg(:), ...
            (m.Const_Machine_R_s(:).*current(:, 1)-m.PE_Ud_meas_avg(:))./m.PE_w_el_rad_avg(:)];
        x = [x; current(keep, :)]; psi = [psi; flux(keep, :)];
        returned = [returned; values(keep, :)]; formula = [formula; source_formula(keep, :)];
        groups = [groups; repmat(g, nnz(keep), 1)];
    end

    %% C: saved merged map. Its shift uses the ACTUAL excitation current.
    projected = [x(:, 1)+k0*x(:, 3), x(:, 2)];
    merged = lookup(map.Id_fine, map.Iq_fine, cat(3, map.Psi_d, map.Psi_q), projected);
    predictions = cat(3, returned, merged);
    own = squeeze(all(isfinite(predictions), 2)); common = all(own, 2);
    own_stats = nan(2, 5); common_stats = nan(2, 5);
    per_group = nan(numel(names), 5, 2);
    for method = 1:2
        own_stats(method, :) = metrics(predictions(:, :, method)-psi, own(:, method));
        common_stats(method, :) = metrics(predictions(:, :, method)-psi, common);
        for g = 1:numel(names)
            per_group(g, :, method) = metrics(predictions(:, :, method)-psi, common & groups == g);
        end
    end

    %% Reconstruct the saved grid from measured points, then vary only raster size.
    F = scatteredInterpolant(projected, psi, 'linear', 'none');
    [D, Q] = meshgrid(map.Id_fine, map.Iq_fine);
    rebuilt = reshape(F([D(:), Q(:)]), [size(D) 2]);
    stored = cat(3, map.Psi_d, map.Psi_q);
    grid_common = all(isfinite(rebuilt) & isfinite(stored), 3);
    grid_mask_difference = nnz(any(isfinite(rebuilt) ~= isfinite(stored), 3));
    grid_difference = reshape(rebuilt-stored, [], 2);
    reconstruction_max = max(abs(grid_difference(grid_common(:), :)));
    node_max = max(abs(F(projected)-psi));
    sizes = [101 201 401]; raster_prediction = nan(size(psi, 1), 2, 3);
    for h = 1:3
        d = linspace(min(map.Id_fine), max(map.Id_fine), sizes(h));
        q = linspace(min(map.Iq_fine), max(map.Iq_fine), sizes(h));
        [D, Q] = meshgrid(d, q);
        values = reshape(F([D(:), Q(:)]), [size(D) 2]);
        raster_prediction(:, :, h) = lookup(d, q, values, projected);
    end
    raster_mask = own(:, 2) & all(isfinite(raster_prediction), [2 3]);
    raster = nan(3, 5);
    for h = 1:3
        raster(h, :) = metrics(raster_prediction(:, :, h)-psi, raster_mask);
    end

    %% Small k sweep: same 101x101 saved grid extent, same points for every k.
    % Reconstruction only: all source samples also train these interpolants.
    factors = sort(unique([12:.5:16 k0]));
    sweep_prediction = nan(size(psi, 1), 2, numel(factors));
    [D, Q] = meshgrid(map.Id_fine, map.Iq_fine);
    for h = 1:numel(factors)
        z = [x(:, 1)+factors(h)*x(:, 3), x(:, 2)];
        G = scatteredInterpolant(z, psi, 'linear', 'none');
        values = reshape(G([D(:), Q(:)]), [size(D) 2]);
        sweep_prediction(:, :, h) = lookup(map.Id_fine, map.Iq_fine, values, z);
    end
    sweep_mask = common & all(isfinite(sweep_prediction), [2 3]); sweep = nan(numel(factors), 5);
    for h = 1:numel(factors)
        sweep(h, :) = metrics(sweep_prediction(:, :, h)-psi, sweep_mask);
    end
    [~, best] = min(sum(sweep(:, 1:2).^2, 2));

    diagnosis{dataset} = struct('file', files{dataset}, 'x', x, 'psi', psi, 'groups', groups, ...
        'excitation', excitation, 'coordinate', coordinate, 'axis_error', axis_error, 'k0', k0, ...
        'prediction', predictions, 'own', own, 'common', common, 'own_stats', own_stats, ...
        'common_stats', common_stats, 'per_group', per_group, ...
        'formula_max_Wb', max(abs(formula-psi)), 'reconstruction_max_Wb', reconstruction_max, ...
        'grid_mask_difference', grid_mask_difference, ...
        'node_max_Wb', node_max, 'unique_projected', size(unique(projected, 'rows'), 1), ...
        'sizes', sizes, 'raster', raster, 'raster_mask', raster_mask, ...
        'factors', factors, 'sweep', sweep, 'sweep_mask', sweep_mask, 'best', best);

    %% One overview per dataset: both components, identical common-point support.
    fig = figure(40+dataset); clf; theme(fig, 'light');
    set(fig, 'Position', [50 30 1100 1300]);
    t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    labels = {'d', 'q'};
    for component = 1:2
        nexttile(component);
        plot(1000*psi(common, component), 1000*returned(common, component), '.', 'MarkerSize', 3); hold on;
        plot(1000*psi(common, component), 1000*merged(common, component), '.', 'MarkerSize', 3);
        limits = 1000*[min(psi(common, component)) max(psi(common, component))];
        plot(limits, limits, 'k--', 'HandleVisibility', 'off'); grid on;
        xlabel('Measured flux in mWb'); ylabel('Reconstructed flux in mWb');
        title(sprintf('Set %g, %s flux', dataset, labels{component}));
        legend('Complete slices', 'Merged 2D', 'Location', 'best');
        errors = 1000*(predictions(:, component, :)-psi(:, component));
        color_limit = max(abs(errors(common, :)), [], 'all');
        for method = 1:2
            nexttile(2*method+component);
            scatter(x(common, 1), x(common, 2), 7, errors(common, 1, method), 'filled');
            clim([-color_limit color_limit]); cb = colorbar; cb.Label.String = 'Flux residual in mWb';
            xlabel('Physical Id in A'); ylabel('Physical Iq in A');
            methods = {'Complete slices', 'Merged 2D'}; title([methods{method} ', ' labels{component} ' flux']); grid on;
        end
        nexttile(6+component);
        plot(excitation, squeeze(per_group(:, component, :)), '.-', 'LineWidth', 1); grid on;
        xlabel('Excitation reference in A'); ylabel('Flux RMSE in mWb');
        title([labels{component} ' flux, common points']); legend('Complete slices', 'Merged 2D', 'Location', 'best');
    end
    exportgraphics(t, fullfile(figures, sprintf('reduction_set%g.png', dataset)), 'Resolution', 180);
    fprintf('\nSet %g: k = %.9f; %g source points; %g unique projected points\n', dataset, k0, size(x, 1), size(unique(projected, 'rows'), 1));
    disp('Rows complete/merged; columns d/q RMSE, d/q maximum in mWb, point count');
    disp('Own support:'); disp(own_stats); disp('Common support:'); disp(common_stats);
    fprintf('Formula, grid rebuild and scattered-node max errors in Wb:\n');
    disp([max(abs(formula-psi)); reconstruction_max; node_max]);
    fprintf('Rebuilt versus saved grid: %g differing finite-mask cells\n', grid_mask_difference);
    disp('Raster 101/201/401, saved extent, common mask:'); disp(raster);
    disp('k sweep: k, d/q RMSE, d/q maximum in mWb, count'); disp([factors' sweep]);
end
save(fullfile(out, 'flux_map_reduction.mat'), 'diagnosis');

fig = figure(43); clf; theme(fig, 'light'); set(fig, 'Position', [50 50 1100 750]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact');
for dataset = 1:2
    a = diagnosis{dataset}; nexttile(dataset);
    plot(a.sizes, a.raster(:, 1:2), '.-', 'LineWidth', 1);
    xlabel('Grid points per axis'); ylabel('Flux RMSE in mWb'); title(sprintf('Set %g, saved extent', dataset));
    legend('d flux', 'q flux'); grid on;
    nexttile(dataset+2); plot(a.factors, a.sweep(:, 1:2), '.-', 'LineWidth', 1); hold on;
    xline(a.k0, 'k--', 'Saved k', 'HandleVisibility', 'off');
    xlabel('k in A/A'); ylabel('Flux RMSE in mWb'); title(sprintf('Set %g, fixed 101 by 101 grid', dataset));
    legend('d flux', 'q flux'); grid on;
end
exportgraphics(t, fullfile(figures, 'reduction_processing.png'), 'Resolution', 180);

function prediction = lookup(d, q, values, z)
    prediction = nan(size(z, 1), 2);
    for component = 1:2
        prediction(:, component) = interp2(d, q, values(:, :, component), z(:, 1), z(:, 2), 'linear', nan);
    end
end

function row = metrics(residual, mask)
    residual = 1000*residual(mask, :);
    row = [sqrt(mean(residual.^2)), max(abs(residual)), nnz(mask)];
end

