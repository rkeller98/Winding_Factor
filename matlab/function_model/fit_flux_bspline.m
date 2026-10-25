function fit_flux_bspline
% Cubic tensor models; sparse weighted QR with a coefficient ridge.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root = setup_matlab_paths('function_model');
out = fullfile(root,'results','flux_fits');
require_analysis_file(fullfile(out,'data.mat'),'prepare_flux_fit_data');
a = load(fullfile(out,'data.mat'));
data = a.data;
models = cell(size(data));
candidates = [];
sizes = [5 5 4;7 7 5;9 9 6];
lambdas = [1e-7 1e-4];

for s = 1:numel(data)
    d = data{s};
    best = inf;

    for j = 1:size(sizes, 1)
        for lambda = lambdas
            n = sizes(j, :);
            m = identify(d, d.train, n, lambda);
            p = eval_flux_model(m, d.x);
            loss = flux_group_loss(p, d.y, d.group, d.validation_support);
            candidates = [candidates;s, n, lambda, 1e6*loss, numeric_bytes(m)];

            if loss < best
                best = loss;
                chosen = n;
                reg = lambda;
            end
        end
    end
    models{s} = identify(d, d.train | d.validation, chosen, reg);
    fprintf('%s spline: %s coefficients/output, lambda %g; validation %.4f mVs\n', ...
        d.name, mat2str(chosen), reg, sqrt(best)*1000);
    save(fullfile(out,'bspline_models.mat'),'models','candidates','sizes','lambdas','-v7');
end
writetable(array2table(candidates,'VariableNames', {'dataset','n_d','n_q','n_third', ...
    'lambda','validation_mVs2','numeric_bytes'}), fullfile(out,'bspline_candidates.csv'));
end

function m = identify(d, use, n, lambda)
x = d.x(use, :);
y = d.y(use, :);
g = d.group(use);
origin = min(x);
scale = max(x)-origin;
z = (x-origin)./scale;
knots = cell(1, 3);

for r = 1:3
    interior = (1:n(r)-4)/(n(r)-3);
    knots{r} = [zeros(1, 4), interior, ones(1, 4)];
end
A = flux_tensor_basis(z, knots);
weight = zeros(size(g));
ug = unique(g);

for j = 1:numel(ug)
    a = g == ug(j);
    weight(a) = 1/sqrt(nnz(a)*numel(ug));
end
C = [A.*weight;sqrt(lambda)*speye(prod(n))]\[y.*weight;zeros(prod(n), 2)];
m = struct('kind','bspline','origin', origin,'scale', scale,'bounds', [min(x);max(x)], ...
    'active', {d.active},'input_units', {d.input_units},'conditions', d.fixed_condition, ...
    'validity','Training bounding box only; unsupported holes require a separate guard', ...
    'knots', {knots},'degree', 3,'coefficients', C,'lambda', lambda,'sizes', n,'training_rows', find(use));
end

function n = numeric_bytes(m)
n = 8*(numel(m.origin)+numel(m.scale)+numel(m.bounds)+sum(cellfun(@numel, m.knots))+1+numel(m.coefficients));
end
