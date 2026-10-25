function loss = flux_group_loss(y, target, group, mask)
% Equal group weight and equal d/q weight; missing predictions have fixed penalty.
r = y-target;
r(~isfinite(r)) = 1;
g = unique(group(mask));
loss = 0;

for j = 1:numel(g)
    a = r(mask & group == g(j), :);
    loss = loss+mean(a.^2,'all')/numel(g);
end
end
