function [study, figures] = investigate_contour_averaging(set_index, group_indices, level_index, export_format)
% Compare contour locations without changing the conventional estimator.
% Default: Set1, a central adjacent group pair, original level 50 of 100.
% group_indices are indices in sorted nominal excitation, level_index is 1:100.
% Example: [study, figures] = investigate_contour_averaging(1, [3 4], 50);
% Optional fourth argument ('png', 'pdf', ...) exports only these figures.
if nargin < 1 || isempty(set_index)
    set_index = 1;
end
if nargin < 2
    group_indices = [];
end
if nargin < 3 || isempty(level_index)
    level_index = 50;
end
if nargin < 4
    export_format = '';
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'matlab', 'common'), fullfile(root, 'matlab', 'k_factor'));
loaded = load(fullfile(root, 'results', 'reduction_data.mat'), 'data');
[currents, flux, groups, excitation] = reduction_inputs(loaded.data{set_index});
flux_span = max(flux(:,1))-min(flux(:,1));
levels = linspace(min(flux(:,1)) + .25*flux_span, ...
    min(flux(:,1)) + .75*flux_span, 100);
level = levels(level_index);
if isempty(group_indices)
    candidates = [];
    for g = 1:numel(excitation)-1
        f1 = flux(groups == g,1);
        f2 = flux(groups == g+1,1);
        if level > max(min(f1),min(f2)) && level < min(max(f1),max(f2))
            candidates(end+1) = g; %#ok<AGROW>
        end
    end
    if isempty(candidates)
        error('No adjacent groups span the selected original level.');
    end
    [~, middle] = min(abs(candidates - numel(excitation)/2));
    group_indices = candidates(middle) + [0 1];
elseif isscalar(group_indices)
    group_indices = group_indices + [0 1];
end
if numel(group_indices) ~= 2 || diff(group_indices) ~= 1
    error('Choose two adjacent sorted excitation-group indices.');
end
study.set = set_index;
study.group_indices = group_indices;
study.nominal_excitation = excitation(group_indices);
study.level_index = level_index;
study.level = level;
study.methods = {'edge mean','length centroid','edge median','edge trim 10%', ...
    'arc sample mean','arc sample median','arc sample trim 10%', ...
    'equal segment mean','thin mesh edge mean','thin mesh length centroid'};
locations = nan(2, numel(study.methods), 3);
departures = nan(2, numel(study.methods));
figures = gobjects(2,1);
figures(1) = figure('Name',sprintf('Set%d contour averaging',set_index), ...
    'Color','w','Position',[70 80 1320 760],'Visible','on');
theme(figures(1), 'light');
layout = tiledlayout(figures(1),2,2,'TileSpacing','compact');
colors = lines(4);
markers = {'o','s','d','^'};
for j = 1:2
    g = group_indices(j);
    c = currents(groups == g,:);
    f = flux(groups == g,1);
    contour = contour_geometry(c,f,level);
    F = scatteredInterpolant(c(:,1:2),f,'linear','none');
    arc = sample_segments(contour,2000,true);
    equal_segment = sample_segments(contour,20,false);
    % Deterministic half-density subsample: ordering is geometric, not acquisition order.
    % Hull points remain; new Delaunay interpolation is explicitly a different field.
    [~, order] = sortrows(c(:,1:2),[2 1]);
    hull = convhull(c(:,1),c(:,2));
    retained = unique([hull; order(1:2:end)]);
    thinned = contour_geometry(c(retained,:),f(retained),level);
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
    locations(j,:,:) = loc;
    departures(j,:) = F(loc(:,1:2)) - level;
    % The thinned centroids have their own flux departure as well as departure in F.
    Ft = scatteredInterpolant(c(retained,1:2),f(retained),'linear','none');
    thinned.own_flux_departure = Ft(loc(9:10,1:2)) - level;
    contour.currents = c;
    contour.flux = f;
    contour.arc_samples = arc;
    contour.equal_segment_samples = equal_segment;
    contour.thinned = thinned;
    % Mirror measured locations in the full convex-hull linear interpolant.
    % This support check does not apply the long-edge filter or test whether
    % the mirrored query lies on the retained contour. Missing mirrors stay missing.
    mirrored_flux = F([c(:,1),-c(:,2)]);
    mirror_valid = isfinite(mirrored_flux);
    contour.symmetry.mirror_supported = mirror_valid;
    contour.symmetry.flux_difference = mirrored_flux-f;
    contour.symmetry.support_fraction = mean(mirror_valid);
    contour.symmetry.rms_flux_difference = sqrt(mean((mirrored_flux(mirror_valid)-f(mirror_valid)).^2));
    contour.symmetry.max_flux_difference = max(abs(mirrored_flux(mirror_valid)-f(mirror_valid)));
    contour.symmetry.contour_iq_range = [min(contour.points(:,2)),max(contour.points(:,2))];
    study.group(j) = contour;
    ax = nexttile(layout,j);
    hold(ax,'on');
    scatter(ax,c(:,1),c(:,2),12,[.78 .78 .78],'filled','DisplayName','measured support');
    if ~isempty(contour.segments)
        plot(ax,[contour.a(:,1),contour.b(:,1)]',[contour.a(:,2),contour.b(:,2)]', ...
            '-','Color',[.15 .15 .15],'HandleVisibility','off');
    end
    plot(ax,contour.points(:,1),contour.points(:,2),'.','Color',[.2 .45 .9], ...
        'DisplayName','actual retained edge intersections');
    if ~isempty(contour.endpoints)
        plot(ax,contour.endpoints(:,1),contour.endpoints(:,2),'rx', ...
            'MarkerSize',10,'LineWidth',1.5,'DisplayName','observed open endpoints');
    end
    for m = 1:4
        plot(ax,loc(m,1),loc(m,2),markers{m},'Color',colors(m,:), ...
            'MarkerFaceColor',colors(m,:),'MarkerSize',8,'DisplayName',study.methods{m});
    end
    axis(ax,'equal');
    grid(ax,'on');
    xlabel(ax,'i_d (A)');
    ylabel(ax,'i_q (A)');
    title(ax,sprintf('Group %d, nominal i_e=%g A; %d intersections, %d open endpoints', ...
        g,excitation(g),size(contour.points,1),size(contour.endpoints,1)));
    legend(ax,'Location','best','FontSize',8);
    ax = nexttile(layout,j+2);
    hold(ax,'on');
    if ~isempty(contour.segments)
        segment_length = contour.lengths;
        histogram(ax,segment_length,20,'FaceColor',[.25 .5 .75]);
        xlabel(ax,'physical contour segment length (A)');
        ylabel(ax,'segment count');
    end
    grid(ax,'on');
    title(ax,{sprintf('Retained length %.2f A; mirrored convex-hull support %.1f%%', ...
        sum(contour.lengths),100*contour.symmetry.support_fraction), ...
        sprintf('RMS flux difference %.4g Vs in the full linear interpolant', ...
        contour.symmetry.rms_flux_difference)});
end
title(layout,sprintf('Original contour level %d/100: psi_d=%.7g Vs (Set%d)', ...
    level_index,level,set_index),'Interpreter','none');
study.locations = locations;
study.flux_departure = departures;
delta = squeeze(locations(2,:,:)-locations(1,:,:));
study.delta_currents = delta;
study.signed_k = -delta(:,1)./delta(:,3);
study.absolute_numerator_k = abs(delta(:,1))./delta(:,3);
study.signed_k_change = study.signed_k-study.signed_k(1);
study.absolute_k_change = study.absolute_numerator_k-study.absolute_numerator_k(1);
study.iq_drift = delta(:,2)./delta(:,3);
% Check the original implementation directly, without recalculating any fits.
[~, reference] = estimate_contour_factor(currents,flux,groups,false,excitation);
study.reference_mean = reshape(reference.means(group_indices,level_index,:),2,3);
study.reference_mean_difference = squeeze(locations(:,1,:))-study.reference_mean;

% Fixed-iq roots are intersections of the retained segments, not independent fits.
% There is no root continuation rule: ambiguous slices remain unresolved.
qlo = max(arrayfun(@(s) min(s.points(:,2)),study.group));
qhi = min(arrayfun(@(s) max(s.points(:,2)),study.group));
q = zeros(0,1);
if qlo <= qhi
    q = linspace(qlo,qhi,201)';
    if qlo <= 0 && qhi >= 0
        q = unique([q;0]);
    end
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
study.fixed_iq.unique_pair = all(counts == 1,2);
study.fixed_iq.signed_k = -fixed_delta(:,1)./fixed_delta(:,3);
study.fixed_iq.absolute_numerator_k = abs(fixed_delta(:,1))./fixed_delta(:,3);
study.fixed_iq.currents = paired;

figures(2) = figure('Name',sprintf('Set%d averaging and fixed-iq comparison',set_index), ...
    'Color','w','Position',[100 110 1320 740],'Visible','on');
theme(figures(2), 'light');
layout = tiledlayout(figures(2),2,2,'TileSpacing','compact');
ax = nexttile(layout);
hold(ax,'on');
plot(ax,1:10,squeeze(locations(:,:,1))','o-');
xticks(ax,1:10);
xticklabels(ax,study.methods);
xtickangle(ax,35);
ylabel(ax,'representative i_d (A)');
grid(ax,'on');
legend(ax,compose('group %d',group_indices));
ax = nexttile(layout);
hold(ax,'on');
plot(ax,1:10,squeeze(locations(:,:,2))','o-');
xticks(ax,1:10);
xticklabels(ax,study.methods);
xtickangle(ax,35);
ylabel(ax,'representative i_q (A)');
grid(ax,'on');
ax = nexttile(layout);
hold(ax,'on');
plot(ax,q,study.fixed_iq.signed_k,'b-','DisplayName','fixed i_q, unique root in each group');
yline(ax,study.signed_k(1),'--','DisplayName','edge mean signed k');
yline(ax,study.signed_k(2),':','DisplayName','length centroid signed k');
xlabel(ax,'shared i_q (A)');
ylabel(ax,'signed k (A/A)');
grid(ax,'on');
legend(ax,'Location','best');
title(ax,sprintf('%d unique pairs / %d slices; %d ambiguous slices', ...
    sum(study.fixed_iq.unique_pair),numel(q),sum(any(counts>1,2))));
ax = nexttile(layout);
hold(ax,'on');
plot(ax,q,counts,'LineWidth',1.2);
xlabel(ax,'shared i_q (A)');
ylabel(ax,'number of distinct contour roots');
grid(ax,'on');
legend(ax,compose('group %d',group_indices),'Location','best');
title(ax,'Zero roots: missing support; multiple roots: no automatic pairing');
title(layout,'Same segments resampled vs thinned mesh (changed interpolant)','Interpreter','none');
drawnow;
fprintf('\nSet%d, groups %d/%d (nominal %.6g/%.6g A), original level %d: %.9g Vs\n', ...
    set_index,group_indices,excitation(group_indices),level_index,level);
for j = 1:2
    fprintf('Group %d: %d points, %d retained intersections, %d segments, %d open endpoints (%d on hull, %d at rejected crossing edges).\n', ...
        group_indices(j),size(study.group(j).currents,1),size(study.group(j).points,1), ...
        numel(study.group(j).lengths),size(study.group(j).endpoints,1), ...
        study.group(j).hull_endpoint_count,study.group(j).filtered_endpoint_count);
    fprintf('Isolated reference edge crossings excluded from segment-length statistics: %d.\n', ...
        study.group(j).isolated_intersection_count);
    fprintf('Mirrored convex-hull interpolant support %.3f%%; RMS/max flux difference %.9g/%.9g Vs; contour iq range [%.6g %.6g] A.\n', ...
        100*study.group(j).symmetry.support_fraction,study.group(j).symmetry.rms_flux_difference, ...
        study.group(j).symmetry.max_flux_difference,study.group(j).symmetry.contour_iq_range);
end
fprintf('Reference mean maximum difference: %.4g A\n',max(abs(study.reference_mean_difference),[],'all'));
fprintf('Method | group1 [id iq ie], flux departure | group2 [id iq ie], flux departure | signed k, abs-numerator k, changes\n');
for m = 1:10
    fprintf('%s | [%.6f %.6f %.6f], %.9g | [%.6f %.6f %.6f], %.9g | %.9f %.9f, %+.9f %+.9f\n', ...
        study.methods{m},locations(1,m,:),departures(1,m),locations(2,m,:),departures(2,m), ...
        study.signed_k(m),study.absolute_numerator_k(m),study.signed_k_change(m),study.absolute_k_change(m));
end
fprintf('Fixed iq: %d/%d unique pairs, %d missing, %d ambiguous slices.\n', ...
    sum(all(counts==1,2)),numel(q),sum(any(counts==0,2)),sum(any(counts>1,2)));
finite_fixed = study.fixed_iq.signed_k(isfinite(study.fixed_iq.signed_k));
if ~isempty(finite_fixed)
    fprintf('Fixed-iq signed k min/max: %.9f / %.9f A/A\n',min(finite_fixed),max(finite_fixed));
end
fprintf('Open endpoints describe the available mesh, not proof of physical contour truncation.\n');
fprintf('Length centroids integrate retained segments only; isolated reference crossings have zero segment length and are omitted.\n');
fprintf('Mirror coverage uses the full convex-hull interpolant, without long-edge filtering; it does not establish retained-contour symmetry.\n');
fprintf('Thinning changes the interpolant; arc/equal-segment resampling only changes averaging weights.\n');
if ~isempty(export_format)
    for i = 1:numel(figures)
        export_analysis_plot(figures(i),fullfile(root,'k_factor_investigation','figures', ...
            sprintf('contour_averaging_set%d_g%d_%d_level%d_%d.png', ...
            set_index,group_indices,level_index,i)),export_format,150);
    end
end
end

function contour = contour_geometry(c,f,level)
% Use exactly the reference mesh, unique edge list and strict long-edge filter.
t = delaunay(c(:,1),c(:,2));
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
degree = accumarray(segments(:),1,[size(nodes,1),1]);
endpoint_nodes = find(degree==1);
hull_edges = freeBoundary(triangulation(t,c(:,1:2)));
on_hull = ismember(edges(kept,:),sort(hull_edges,2),'rows');
contour.points = points;
contour.edges = edges(kept,:);
contour.nodes = nodes;
contour.segments = segments;
contour.a = a;
contour.b = b;
contour.lengths = lengths;
% Exact line integral over retained segments; isolated crossings contribute no length.
contour.length_centroid = sum(lengths.*((a+b)/2),1)/sum(lengths);
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
