function fit_flux_rbf
% Configurable compact Gaussian models, fixed bases and linear coefficients.
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'common'));
root = setup_matlab_paths('function_model');
out = fullfile(root,'results','flux_fits');
require_analysis_file(fullfile(out,'data.mat'),'prepare_flux_fit_data');
a = load(fullfile(out,'data.mat'));
data = a.data;
models = cell(size(data));
candidates = [];
center_counts = [32 64 128];
width_factors = [.75 1.5 3];
lambdas = [1e-7 1e-4];

for s = 1:numel(data)
    d = data{s};
    best = inf;

    for count = center_counts
        for width = width_factors
            for lambda = lambdas
                m = identify(d, d.train, count, width, lambda);
                p = eval_flux_model(m, d.x);
                loss = flux_group_loss(p, d.y, d.group, d.validation_support);
                candidates = [candidates;s, count, width, lambda, 1e6*loss, model_bytes(m)];

                if loss < best
                    best = loss;
                    chosen = [count, width, lambda];
                end
            end
        end
    end
    models{s} = identify(d, d.train | d.validation, chosen(1), chosen(2), chosen(3));
    fprintf('%s RBF: %d centers, width factor %.2g, lambda %g; validation %.4f mVs\n', ...
        d.name, chosen, sqrt(best)*1000);
    save(fullfile(out,'rbf_models.mat'),'models','candidates','center_counts','width_factors','lambdas','-v7');
end
writetable(array2table(candidates,'VariableNames', {'dataset','centers','width_factor', ...
    'lambda','validation_mVs2','numeric_bytes'}), fullfile(out,'rbf_candidates.csv'));
end

function m = identify(d, use, count, width, lambda)
x = d.x(use, :);
y = d.y(use, :);
g = d.group(use);
origin = min(x);
scale = max(x)-origin;
z = (x-origin)./scale;
count = min(count, size(z, 1));
centers = zeros(count, size(z, 2));
[~, ix] = min(sum((z-mean(z)).^2, 2));
distance = inf(size(z, 1), 1);

for j = 1:count
    centers(j, :) = z(ix, :);
    distance = min(distance, sum((z-z(ix, :)).^2, 2));
    [~, ix] = max(distance);
end
dist = sum((reshape(centers, count, 1, [])-reshape(centers, 1, count, [])).^2, 3);
dist(1:count+1:end) = inf;
sigma = width*median(sqrt(min(dist, [], 2)));
phi = exp(-sum((reshape(z, size(z, 1), 1, [])-reshape(centers, 1, count, [])).^2, 3)/sigma^2);
A = [ones(size(z, 1), 1), z, phi];
weight = zeros(size(g));
ug = unique(g);

for j = 1:numel(ug)
    a = g == ug(j);
    weight(a) = 1/sqrt(nnz(a)*numel(ug));
end
penalty = diag([zeros(1, size(z, 2)+1), ones(1, count)]);
C = [A.*weight;sqrt(lambda)*penalty]\[y.*weight;zeros(size(A, 2), 2)];
m = struct('kind','rbf','origin', origin,'scale', scale,'bounds', [min(x);max(x)], ...
    'active', {d.active},'input_units', {d.input_units},'conditions', d.fixed_condition, ...
    'validity','Training bounding box only; unsupported holes require a separate guard', ...
    'centers', centers,'width', sigma,'coefficients', C,'lambda', lambda, ...
    'width_factor', width,'training_rows', find(use));
end

function n = model_bytes(m)
n = 8*(numel(m.origin)+numel(m.scale)+numel(m.bounds)+numel(m.centers)+1+numel(m.coefficients));
end
