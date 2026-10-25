function study = contour_pair_statistics(currents, flux, groups, group_indices, level, iq_queries, meshes)
% Pure numerical contour comparison; group_indices index sorted unique(groups).
% Optional iq_queries are shared physical q currents, never extrapolated.
% Optional meshes is a two-element cell of per-group triangle connectivity.
% Empty meshes use the unchanged reference Delaunay convention.
if nargin < 6
    iq_queries = [];
end
if nargin < 7
    meshes = {[], []};
end
study.group_indices = group_indices;
study.level = level;
study.locations = nan(2,10,3);
study.flux_departure = nan(2,10);
study.mean_iq_squared = nan(2,10);
group_ids = unique(groups);
for j = 1:2
    rows = groups == group_ids(group_indices(j));
    item = contour_group_statistics(currents(rows,:),flux(rows,1),level,meshes{j});
    study.methods = item.methods;
    study.locations(j,:,:) = item.locations;
    study.flux_departure(j,:) = item.flux_departure;
    study.mean_iq_squared(j,:) = item.mean_iq_squared;
    study.group(j) = item.geometry;
end
delta = squeeze(study.locations(2,:,:)-study.locations(1,:,:));
study.delta_currents = delta;
study.signed_k = -delta(:,1)./delta(:,3);
study.absolute_numerator_k = abs(delta(:,1))./delta(:,3);
study.signed_k_change = study.signed_k-study.signed_k(1);
study.absolute_k_change = study.absolute_numerator_k-study.absolute_numerator_k(1);
study.iq_drift = delta(:,2)./delta(:,3);
study.fixed_iq = contour_fixed_iq(study.group,iq_queries);
end
