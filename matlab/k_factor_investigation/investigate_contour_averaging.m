function [study, figures] = investigate_contour_averaging(set_index, group_indices, level_index, export_format)
% Compare contour locations without changing the conventional estimator.
% Default: Set1, a central adjacent group pair, original level 50 of 100.
% group_indices are indices in sorted nominal excitation, level_index is 1:100.
% Example: [study, figures] = investigate_contour_averaging(1, [3 4], 50);
% Also compares original levels 25, 50 and 75 for the same group pair.
% Optional fourth argument ('png', 'pdf', ...) exports into the local results folder.
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
study = contour_pair_statistics(currents,flux,groups,group_indices,level);
study.set = set_index;
study.nominal_excitation = excitation(group_indices);
study.level_index = level_index;
locations = study.locations;
departures = study.flux_departure;
figures = gobjects(2,1);
figures(1) = figure('Name',sprintf('Set%d contour averaging',set_index), ...
    'Color','w','Position',[70 80 1320 760],'Visible','on');
theme(figures(1), 'light');
layout = tiledlayout(figures(1),2,2,'TileSpacing','compact');
colors = lines(4);
markers = {'o','s','d','^'};
for j = 1:2
    g = group_indices(j);
    contour = study.group(j);
    c = contour.currents;
    loc = squeeze(locations(j,:,:));
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
% Keep the original estimator as an independent numerical reference.
[~,reference] = estimate_contour_factor(currents,flux,groups,false,excitation);
study.reference_mean = reshape(reference.means(group_indices,level_index,:),2,3);
study.reference_mean_difference = squeeze(locations(:,1,:))-study.reference_mean;
q = study.fixed_iq.q;
counts = study.fixed_iq.root_counts;
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
% Representative level spread, without refitting Method B or any flux model.
comparison_indices = [25 50 75];
comparison_k = nan(3,numel(study.methods));
fixed_ranges = nan(3,3);
coverage = nan(3,1);
q_ranges = nan(3,2);
for i = 1:numel(comparison_indices)
    index = comparison_indices(i);
    comparison = contour_pair_statistics(currents,flux,groups,group_indices,levels(index));
    comparison_k(i,:) = comparison.signed_k';
    fixed = comparison.fixed_iq.signed_k;
    fixed = fixed(isfinite(fixed));
    if ~isempty(fixed)
        fixed_ranges(i,:) = [min(fixed),median(fixed),max(fixed)];
    end
    coverage(i) = mean(comparison.fixed_iq.unique_pair);
    if ~isempty(comparison.fixed_iq.q)
        q_ranges(i,:) = [min(comparison.fixed_iq.q),max(comparison.fixed_iq.q)];
    end
    fprintf('Original level %d, %.7g Vs: signed edge/length/median/trim k %.6f %.6f %.6f %.6f; fixed-q min/median/max %.6f %.6f %.6f; unique fraction %.1f%% within queried overlap [%.6g %.6g] A.\n', ...
        index,levels(index),comparison_k(i,1:4),fixed_ranges(i,:),100*coverage(i),q_ranges(i,:));
end
study.level_comparison.level_indices = comparison_indices;
study.level_comparison.levels = levels(comparison_indices);
study.level_comparison.methods = study.methods;
study.level_comparison.signed_k = comparison_k;
study.level_comparison.fixed_iq_min_median_max = fixed_ranges;
study.level_comparison.unique_coverage = coverage;
study.level_comparison.queried_iq_ranges = q_ranges;
fprintf('All methods on the three original levels (A/A):\n');
disp(array2table(comparison_k,'VariableNames',matlab.lang.makeValidName(study.methods), ...
    'RowNames',compose('level_%d',comparison_indices)));
figures(3) = figure('Name',sprintf('Set%d contour methods across levels',set_index), ...
    'Color','w','Position',[130 140 1050 600],'Visible','on');
theme(figures(3),'light');
ax = axes(figures(3));
hold(ax,'on');
plot(ax,levels(comparison_indices),comparison_k(:,1:4),'o-');
errorbar(ax,levels(comparison_indices),fixed_ranges(:,2), ...
    fixed_ranges(:,2)-fixed_ranges(:,1),fixed_ranges(:,3)-fixed_ranges(:,2),'ks-', ...
    'MarkerFaceColor','k');
xlabel(ax,'original contour level psi_d (Vs)');
ylabel(ax,'signed k (A/A)');
grid(ax,'on');
legend(ax,[study.methods(1:4),{'fixed-q median and min/max over q'}],'Location','best');
title(ax,{'Three original levels; measured i_e differences used throughout', ...
    'Direct fixed-q contour roots differ from existing Method B prediction optimization'});
fprintf('Fixed-q factors use interpolated measured i_e, not nominal group-reference differences.\n');
fprintf('Direct fixed-q contour matching is distinct from existing Method B fixed-q prediction optimization.\n');
drawnow;
if ~isempty(export_format)
    for i = 1:numel(figures)
        export_analysis_plot(figures(i),fullfile(root,'matlab','k_factor_investigation','results', ...
            sprintf('contour_averaging_set%d_g%d_%d_level%d_%d.png', ...
            set_index,group_indices,level_index,i)),export_format,150);
    end
end
end
