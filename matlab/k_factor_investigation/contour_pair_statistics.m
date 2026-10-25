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
study.methods = {'edge mean','length centroid','edge median','edge trim 10%', ...
    'arc sample mean','arc sample median','arc sample trim 10%', ...
    'equal segment mean','thin mesh edge mean','thin mesh length centroid'};
study.locations = nan(2,numel(study.methods),3);
study.flux_departure = nan(2,numel(study.methods));
% Robust componentwise locations do not share one set of averaging weights.
% Their second moments remain NaN rather than suggesting the mean identity applies.
study.mean_iq_squared = nan(2,numel(study.methods));
group_ids = unique(groups);
for j = 1:2
    rows = groups == group_ids(group_indices(j));
    c = currents(rows,:);
    f = flux(rows,1);
    contour = contour_geometry(c,f,level,meshes{j});
    F = scatteredInterpolant(c(:,1:2),f,'linear','none');
    arc = sample_segments(contour,2000,true);
    equal_segment = sample_segments(contour,20,false);
    % Hull-preserving thinning changes the interpolant, not just sample weights.
    [~,order] = sortrows(c(:,1:2),[2 1]);
    hull = convhull(c(:,1),c(:,2));
    retained = unique([hull;order(1:2:end)]);
    thinned = contour_geometry(c(retained,:),f(retained),level,[]);
    loc = [
        mean(contour.points,1)
        contour.length_centroid
        median(contour.points,1)
        trimmed_location(contour.points)
        mean(arc,1)
        median(arc,1)
        trimmed_location(arc)
        mean(equal_segment,1)
        mean(thinned.points,1)
        thinned.length_centroid
    ];
    study.locations(j,:,:) = loc;
    study.flux_departure(j,:) = F(loc(:,1:2))-level;
    if ~isempty(meshes{j})
        mesh = triangulation(meshes{j},c(:,1:2));
        [triangle,barycentric] = pointLocation(mesh,loc(:,1:2));
        valid = isfinite(triangle);
        values = nan(size(loc,1),1);
        vertices = meshes{j}(triangle(valid),:);
        values(valid) = sum(barycentric(valid,:).*f(vertices),2);
        study.flux_departure(j,:) = values-level;
    end
    study.mean_iq_squared(j,[1 2 5 8 9 10]) = [
        mean(contour.points(:,2).^2)
        contour.length_mean_iq_squared
        mean(arc(:,2).^2)
        mean(equal_segment(:,2).^2)
        mean(thinned.points(:,2).^2)
        thinned.length_mean_iq_squared
    ];
    Ft = scatteredInterpolant(c(retained,1:2),f(retained),'linear','none');
    thinned.own_flux_departure = Ft(loc(9:10,1:2))-level;
    contour.currents = c;
    contour.flux = f;
    contour.delaunay_crossing_residual = F(contour.points(:,1:2))-level;
    contour.arc_samples = arc;
    contour.equal_segment_samples = equal_segment;
    contour.thinned = thinned;
    % Full convex-hull interpolant support; the long-edge filter is absent here.
    mirrored_flux = F([c(:,1),-c(:,2)]);
    mirror_valid = isfinite(mirrored_flux);
    contour.symmetry.mirror_supported = mirror_valid;
    contour.symmetry.flux_difference = mirrored_flux-f;
    contour.symmetry.support_fraction = mean(mirror_valid);
    contour.symmetry.rms_flux_difference = sqrt(mean((mirrored_flux(mirror_valid)-f(mirror_valid)).^2));
    contour.symmetry.max_flux_difference = max(abs(mirrored_flux(mirror_valid)-f(mirror_valid)));
    contour.symmetry.contour_iq_range = [min(contour.points(:,2)),max(contour.points(:,2))];
    study.group(j) = contour;
end
delta = squeeze(study.locations(2,:,:)-study.locations(1,:,:));
study.delta_currents = delta;
study.signed_k = -delta(:,1)./delta(:,3);
study.absolute_numerator_k = abs(delta(:,1))./delta(:,3);
study.signed_k_change = study.signed_k-study.signed_k(1);
study.absolute_k_change = study.absolute_numerator_k-study.absolute_numerator_k(1);
study.iq_drift = delta(:,2)./delta(:,3);
% No branch continuation: retain all roots and use only unique-root pairs.
if isempty(iq_queries)
    q = zeros(0,1);
    if ~any(arrayfun(@(s) isempty(s.points),study.group))
        qlo = max(arrayfun(@(s) min(s.points(:,2)),study.group));
        qhi = min(arrayfun(@(s) max(s.points(:,2)),study.group));
        if qlo <= qhi
            q = linspace(qlo,qhi,201)';
            if qlo <= 0 && qhi >= 0
                q = unique([q;0]);
            end
        end
    end
else
    q = iq_queries(:);
end
roots = cell(numel(q),2);
counts = zeros(numel(q),2);
paired = nan(numel(q),2,3);
for i = 1:numel(q)
    for j = 1:2
        roots{i,j} = fixed_iq_roots(study.group(j),q(i));
        counts(i,j) = size(roots{i,j},1);
        if counts(i,j) == 1
            paired(i,j,:) = roots{i,j};
        end
    end
end
fixed_delta = reshape(paired(:,2,:)-paired(:,1,:),numel(q),3);
study.fixed_iq.q = q;
study.fixed_iq.roots = roots;
study.fixed_iq.root_counts = counts;
study.fixed_iq.unique_pair = all(counts==1,2);
study.fixed_iq.signed_k = -fixed_delta(:,1)./fixed_delta(:,3);
study.fixed_iq.absolute_numerator_k = abs(fixed_delta(:,1))./fixed_delta(:,3);
study.fixed_iq.currents = paired;
end
function contour = contour_geometry(c,f,level,t)
% Use exactly the reference mesh, unique edge list and strict long-edge filter.
if isempty(t)
    t = delaunay(c(:,1),c(:,2));
end
[edges,~,edge_index] = unique(sort([t(:,[1 2]);t(:,[1 3]);t(:,[2 3])],2),'rows');
triangle_edges = reshape(edge_index,size(t,1),3);
length_squared = sum((c(edges(:,1),1:2)-c(edges(:,2),1:2)).^2,2);
short = length_squared < 9*mean(length_squared);
w = (level-f(edges(:,1)))./(f(edges(:,2))-f(edges(:,1)));
crossing = w>=0 & w<=1;
kept = find(crossing & short);
points = c(edges(kept,1),:) + w(kept).*(c(edges(kept,2),:)-c(edges(kept,1),:));
% Mapping edge IDs to distinct geometric nodes handles shared contour vertices.
[nodes,~,node_ids] = unique(points(:,1:2),'rows');
edge_node = zeros(size(edges,1),1);
edge_node(kept) = node_ids;
segments = zeros(0,2);
a = zeros(0,3);
b = zeros(0,3);
filtered_nodes = [];
for k = 1:size(t,1)
    e = triangle_edges(k,:);
    n = unique(edge_node(e));
    n(n==0) = [];
    if numel(n)==2
        p1 = points(find(node_ids==n(1),1),:);
        p2 = points(find(node_ids==n(2),1),:);
        segments(end+1,:) = n(:)'; %#ok<AGROW>
        a(end+1,:) = p1; %#ok<AGROW>
        b(end+1,:) = p2; %#ok<AGROW>
    elseif any(crossing(e) & ~short(e))
        filtered_nodes = [filtered_nodes;n(:)]; %#ok<AGROW>
    end
end
% An iso-edge shared by triangles can otherwise appear twice.
[segments,ia] = unique(sort(segments,2),'rows');
a = a(ia,:);
b = b(ia,:);
lengths = sqrt(sum((b(:,1:2)-a(:,1:2)).^2,2));
degree = zeros(size(nodes,1),1);
if ~isempty(segments)
    degree = accumarray(segments(:),1,[size(nodes,1),1]);
end
endpoint_nodes = find(degree==1);
hull_edges = freeBoundary(triangulation(t,c(:,1:2)));
on_hull = ismember(edges(kept,:),sort(hull_edges,2),'rows');
contour.points = points;
contour.triangles = t;
contour.edges = edges(kept,:);
contour.nodes = nodes;
contour.segments = segments;
contour.a = a;
contour.b = b;
contour.lengths = lengths;
% Exact line integral over retained segments; isolated crossings contribute no length.
contour.length_centroid = sum(lengths.*((a+b)/2),1)/sum(lengths);
contour.length_mean_iq_squared = sum(lengths.*(a(:,2).^2+a(:,2).*b(:,2)+b(:,2).^2)/3)/sum(lengths);
contour.endpoints = nodes(endpoint_nodes,:);
contour.hull_endpoint_count = sum(ismember(endpoint_nodes,node_ids(on_hull)));
contour.filtered_endpoint_count = sum(ismember(endpoint_nodes,filtered_nodes));
contour.isolated_intersection_count = sum(degree(node_ids)==0);
contour.rejected_crossing_edge_count = sum(crossing & ~short);
end

function location = trimmed_location(points)
% Componentwise trim: remove floor(10% N) values independently at each tail.
ordered = sort(points,1);
n = size(points,1);
trim = floor(.1*n);
location = mean(ordered(1+trim:n-trim,:),1);
end

function points = sample_segments(contour,n,by_length)
% Midpoint quadrature: arc samples have equal physical-length weights.
% Equal-segment samples give every segment n samples irrespective of length.
points = zeros(0,3);
if isempty(contour.lengths)
    return;
end
if by_length
    cumulative = [0;cumsum(contour.lengths)];
    distance = ((.5:n-.5)'/n)*cumulative(end);
    [~,segment] = histc(distance,cumulative);
    u = (distance-cumulative(segment))./contour.lengths(segment);
    points = contour.a(segment,:) + u.*(contour.b(segment,:)-contour.a(segment,:));
else
    u = ((.5:n-.5)'/n);
    for i = 1:numel(contour.lengths)
        points = [points;contour.a(i,:)+u.*(contour.b(i,:)-contour.a(i,:))]; %#ok<AGROW>
    end
end
end

function roots = fixed_iq_roots(contour,q)
a = contour.a;
b = contour.b;
u = (q-a(:,2))./(b(:,2)-a(:,2));
crossing = u>=0 & u<=1;
roots = a(crossing,:) + u(crossing).*(b(crossing,:)-a(crossing,:));
% A horizontal segment represents a continuum of roots: retain both ends,
% so the slice is explicitly ambiguous and never enters an automatic pair.
horizontal = a(:,2)==q & b(:,2)==q;
roots = [roots;a(horizontal,:);b(horizontal,:)];
if ~isempty(roots)
    roots(:,2) = q;
    roots = uniquetol(roots,1e-10,'ByRows',true,'DataScale',1);
    roots = sortrows(roots,1);
end
end
