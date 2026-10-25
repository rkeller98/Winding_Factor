function map = fit_reduction_map(x, y, model, p, current_scale, n, normalize)
% Linear scattered field -> stored grid. Reject long simplices, no extrapolation.
z = reduction_coordinates(x, model, p, current_scale);
origin = min(z);
sc = max(z)-origin;
if ~normalize
    sc = ones(size(sc));
end
z = (z-origin)./sc;
[z, ~, ic] = unique(z,'rows');
v = zeros(size(z, 1), 2);

for c = 1:2
    v(:, c) = accumarray(ic, y(:, c), [], @mean);
end
dt = delaunayTriangulation(z);
t = dt.ConnectivityList;
comb = nchoosek(1:size(t, 2), 2);
length2 = zeros(size(t, 1), size(comb, 1));

for j = 1:size(comb, 1)
    length2(:, j) = sum((z(t(:, comb(j, 1)), :)-z(t(:, comb(j, 2)), :)).^2, 2);
end
edge = edges(dt);
limit = 9*mean(sum((z(edge(:, 1), :)-z(edge(:, 2), :)).^2, 2));
good = all(length2 < limit, 2);
dim = size(z, 2);
axes = cell(1, dim);

for j = 1:dim
    count = n;
    if j == 3
        count = max(5, numel(unique(round(x(:, 3)))));
        count = min(count, 21);
    end
    axes{j} = linspace(min(z(:, j)), max(z(:, j)), count);
end
grid = cell(1, dim);
[grid{:}] = ndgrid(axes{:});
query = [];
for j = 1:dim
    query = [query, grid{j}(:)];
end
[which, bary] = pointLocation(dt, query);
valid = isfinite(which);
valid(valid) = good(which(valid));
value = nan(size(query, 1), 2);

for c = 1:2
    yy = v(:, c);
    value(valid, c) = sum(bary(valid, :).*yy(t(which(valid), :)), 2);
end
values = cell(1, 2);
for c = 1:2
    values{c} = reshape(value(:, c), size(grid{1}));
end
map = struct('model', model,'p', p,'current_scale', current_scale,'origin', origin, ...
    'scale', sc,'axes', {axes},'values', {values},'nodes', numel(grid{1}), ...
    'valid_nodes', nnz(valid),'source_points', size(z, 1),'long_simplex_limit', limit);
end
