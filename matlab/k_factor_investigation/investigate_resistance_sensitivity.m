function [study, figures] = investigate_resistance_sensitivity(set_index, offsets, group_indices, export_format, iq_query)
% Add constant Rs offsets to ORIGINAL pointwise resistance; baseline stays intact.
% Offsets are fractions of median original Rs, default [-0.1 0 0.1].
% Optional fifth argument fixes the shared physical iq; default is nonzero.
if nargin < 1 || isempty(set_index)
    set_index = 1;
end
if nargin < 2 || isempty(offsets)
    offsets = [-.1 0 .1];
end
if nargin < 3
    group_indices = [];
end
if nargin < 4
    export_format = '';
end
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(fullfile(root,'matlab','common'),fullfile(root,'matlab','k_factor'));
loaded = load(fullfile(root,'results','reduction_data.mat'),'data');
data = loaded.data{set_index};
study = struct('set',set_index,'skipped',true);
figures = gobjects(0);
needed = {'ud_meas','uq_meas','omega_el','stator_resistance','id_meas','iq_meas'};
for j = 1:numel(needed)
    name = needed{j};
    if ~isfield(data,name) || isempty(data.(name)) || any(~isfinite(data.(name)))
        fprintf('Set%d Rs sensitivity skipped: missing/nonfinite %s channel.\n',set_index,name);
        return;
    end
end
if any(abs(data.omega_el) <= 1)
    fprintf('Set%d Rs sensitivity skipped: stationary flux needs |omega_el| > 1 rad/s at every canonical point.\n',set_index);
    return;
end
[x,y,groups,excitation] = reduction_inputs(data);
R = data.stator_resistance;
w = data.omega_el;
recomputed = [(data.uq_meas-R.*x(:,2))./w,(R.*x(:,1)-data.ud_meas)./w];
change_per_ohm = [-x(:,2)./w,x(:,1)./w];
baseline_error = recomputed-y;
study.zero_offset_max_flux_difference = max(abs(baseline_error),[],'all');
fprintf('Set%d: zero-offset voltage recomputation minus original flux: max %.9g Vs.\n', ...
    set_index,study.zero_offset_max_flux_difference);
fprintf('Original pointwise Rs range %.9g--%.9g ohm, median %.9g ohm.\n',min(R),max(R),median(R));
flux_span = max(y(:,1))-min(y(:,1));
levels = linspace(min(y(:,1))+.25*flux_span,min(y(:,1))+.75*flux_span,100);
level = levels(50);
if isempty(group_indices)
    candidates = [];
    for g = 1:numel(excitation)-1
        f1 = y(groups==g,1);
        f2 = y(groups==g+1,1);
        if level > max(min(f1),min(f2)) && level < min(max(f1),max(f2))
            candidates(end+1) = g; %#ok<AGROW>
        end
    end
    if isempty(candidates)
        fprintf('Set%d: no adjacent pair spans original level 50; Rs study skipped.\n',set_index);
        return;
    end
    [~,j] = min(abs(candidates-numel(excitation)/2));
    group_indices = candidates(j)+[0 1];
end
original_pair = contour_pair_statistics(x,y,groups,group_indices,level);
eligible = find(original_pair.fixed_iq.unique_pair & isfinite(original_pair.fixed_iq.signed_k));
if nargin < 5 || isempty(iq_query)
    if isempty(eligible)
        iq_query = 0;
    else
        available_q = original_pair.fixed_iq.q(eligible);
        target = .5*max(available_q);
        [~,j] = min(abs(available_q-target));
        iq_query = available_q(j);
    end
end
original_fixed = contour_pair_statistics(x,y,groups,group_indices,level,iq_query);
baseline_fixed_valid = original_fixed.fixed_iq.unique_pair(1);
if ~baseline_fixed_valid
    fprintf('No unique original root pair at the chosen iq; fixed-iq comparison remains NaN.\n');
end
saved = load(fullfile(root,'results','winding_validation.mat'),'study');
L0 = saved.study{set_index}.slopes{2}.L;
% Fit only the flux corrections. Affine fitting is linear in its output data.
% Use exactly the saved 80-neighbour coordinates, weights and QR screen.
scale = max(x)-min(x);
z = x./scale;
correction = [baseline_error change_per_ohm];
dL = nan(size(x,1),4,3);
for i = 1:size(x,1)
    [distance,ix] = mink(sum((z-z(i,:)).^2,2),80);
    if numel(unique(groups(ix))) < 2
        continue;
    end
    weight = sqrt(exp(-distance/max(distance)));
    A = [ones(80,1),z(ix,:)-z(i,:)].*weight;
    [~,r] = qr(A,0);
    if rcond(r) < 1e-6
        continue;
    end
    coefficients = A\(correction(ix,:).*weight);
    dL(i,:,:) = permute(coefficients(2:4,:)./scale',[3 2 1]);
end
offsets = unique([offsets(:);0])';
delta_R = offsets*median(R);
study.skipped = false;
study.offset_fractions = offsets;
study.delta_R = delta_R;
study.group_indices = group_indices;
study.nominal_excitation = excitation(group_indices);
study.fixed_level = level;
study.fixed_iq = iq_query;
study.original_pair = original_pair;
study.original_fixed = original_fixed;
study.summary = nan(numel(offsets),12);
fprintf('Pair-only comparison holds original lambda %.9g Vs and iq %.9g A fixed (groups %d/%d).\n',level,iq_query,group_indices);
fprintf('Global contour k retains the 100-level middle-range rule, applied to each recomputed flux map.\n');
for j = 1:numel(offsets)
    varied_R = R+delta_R(j);
    flux = [(data.uq_meas-varied_R.*x(:,2))./w,(varied_R.*x(:,1)-data.ud_meas)./w];
    [k,contour] = estimate_contour_factor(x,flux,groups,false,excitation);
    pair = contour_pair_statistics(x,flux,groups,group_indices,level,iq_query);
    varied = data;
    varied.psi_d = flux(:,1);
    varied.psi_q = flux(:,2);
    symmetry = investigate_flux_symmetry(set_index,group_indices,'',varied,false);
    L = L0+dL(:,1:2,:)+delta_R(j)*dL(:,3:4,:);
    a = L(:,:,1);
    b = L(:,:,3);
    an = vecnorm(a,2,2);
    valid = all(isfinite([a b]),2) & an > .1*median(an,'omitnan');
    klocal = sum(a.*b,2)./sum(a.^2,2);
    klocal(~valid) = nan;
    ratios = b./a;
    component_valid = false(size(ratios));
    for c = 1:2
        component_valid(:,c) = isfinite(ratios(:,c)) & abs(a(:,c)) >= .1*median(abs(a(:,c)),'omitnan');
        ratios(~component_valid(:,c),c) = nan;
    end
    kg = sum(a(valid,:).*b(valid,:),'all')/sum(a(valid,:).^2,'all');
    reciprocity = L(:,1,2)-L(:,2,1);
    item = struct('flux',flux,'k',k,'contour',contour,'pair',pair,'symmetry',symmetry, ...
        'L',L,'klocal',klocal,'valid',valid,'component_factors',ratios, ...
        'component_valid',component_valid,'global_derivative_k',kg,'reciprocity',reciprocity);
    study.cases(j) = item;
    fixed_k = pair.fixed_iq.signed_k(1);
    if ~baseline_fixed_valid
        fixed_k = nan;
    end
    sym_d = sqrt(mean(symmetry.summary(:,4).^2,'omitnan'));
    sym_q = sqrt(mean(symmetry.summary(:,5).^2,'omitnan'));
    study.summary(j,:) = [delta_R(j),k,pair.signed_k(1),fixed_k,pair.iq_drift(1), ...
        sym_d,sym_q,median(klocal,'omitnan'),kg,median(reciprocity,'omitnan'), ...
        median(ratios(:,1),'omitnan'),median(ratios(:,2),'omitnan')];
end
common = true(size(groups));
component_common = true(numel(groups),2);
for j = 1:numel(offsets)
    common = common & study.cases(j).valid;
    component_common = component_common & study.cases(j).component_valid;
end
study.common_scalar = common;
study.common_components = component_common;
for j = 1:numel(offsets)
    item = study.cases(j);
    fprintf('delta Rs=%+.7g ohm: global contour k %.8g; fixed-level edge/fixed-q %.8g/%.8g; median local LS %.8g (%d common points).\n', ...
        delta_R(j),item.k,item.pair.signed_k(1),study.summary(j,4),median(item.klocal(common),'omitnan'),nnz(common));
    fprintf('  common-component median k(d)/k(q) %.8g/%.8g; median Ldq-Lqd %.8g mH.\n', ...
        median(item.component_factors(component_common(:,1),1),'omitnan'), ...
        median(item.component_factors(component_common(:,2),2),'omitnan'),1000*median(item.reciprocity,'omitnan'));
end
% Constant-speed prediction is a reference only; varying omega adds derivatives.
reciprocity_slope = dL(:,3,2)-dL(:,4,1);
study.reciprocity_change_per_ohm = reciprocity_slope;
study.constant_speed_reference_per_ohm = -2./w;
fprintf('Derivative reciprocity slope median %.9g H/ohm; -2/omega median %.9g H/ohm.\n', ...
    median(reciprocity_slope,'omitnan'),median(-2./w));
fprintf('Measured omega varies from %.9g to %.9g rad/s; its variation and local-affine estimation add terms to -2 deltaRs/omega.\n',min(w),max(w));
fprintf('This is Rs sensitivity, not Rs identification. Fixed-q is one fixed-level pair, not global Method B.\n');
fprintf('At constant omega, fixed-q comparisons shift the effective original flux level by deltaRs*iq/omega; invariance is not guaranteed.\n');
fprintf('iq=0 is a special low-sensitivity case; varying measured omega and linear interpolation can also matter.\n');
study.summary_columns = {'delta_R_ohm','global_contour_k','fixed_level_edge_signed_k','fixed_q_signed_k', ...
    'centroid_iq_drift','d_mirror_residual_rms_Vs','q_mirror_residual_rms_Vs', ...
    'median_local_LS','global_derivative_k','median_reciprocity_H','median_component_d','median_component_q'};
[fig,plot_tab] = investigation_tab(set_index,'Statorwiderstand','Sensitivität');
selection = uitabgroup(plot_tab,'Position',[0 .94 1 .06]);
for g = 1:numel(excitation)-1
    tab = uitab(selection,'Title',sprintf('%g -> %g A',excitation(g:g+1)),'UserData',g);
    if g==group_indices(1)
        selection.SelectedTab = tab;
    end
end
selection.SelectionChangedFcn = @(~,event) investigate_resistance_sensitivity( ...
    set_index,offsets,event.NewValue.UserData+[0 1],'',[]);
panel = uipanel(plot_tab,'BorderType','none','Position',[0 0 1 .94]);
layout = tiledlayout(panel,2,3,'TileSpacing','compact');
ax = nexttile(layout);
plot(ax,delta_R,study.summary(:,2),'o-');
ylabel(ax,'global contour k (A/A)');
ax = nexttile(layout);
plot(ax,delta_R,study.summary(:,3:4),'o-');
ylabel(ax,'signed pair k (A/A)');
legend(ax,'edge mean','fixed iq','Location','best');
ax = nexttile(layout);
hold(ax,'on');
for g = 1:2
    locations = arrayfun(@(a) a.pair.locations(g,1,2),study.cases);
    plot(ax,delta_R,locations,'o-','DisplayName',sprintf('%g A nominal',excitation(group_indices(g))));
end
ylabel(ax,'edge centroid i_q (A)');
legend(ax,'Location','best');
ax = nexttile(layout);
plot(ax,delta_R,1000*study.summary(:,6:7),'o-');
ylabel(ax,'mirror residual RMS (mVs)');
legend(ax,'d difference','q sum','Location','best');
ax = nexttile(layout);
hold(ax,'on');
for j = 1:numel(offsets)
    plot(ax,delta_R(j),median(study.cases(j).klocal(common),'omitnan'),'bo');
end
ylabel(ax,'local LS median, common support (A/A)');
ax = nexttile(layout);
plot(ax,delta_R,1000*study.summary(:,10),'o-');
ylabel(ax,'median L_{dq}-L_{qd} (mH)');
for ax = findall(plot_tab,'Type','axes')'
    xlabel(ax,'constant delta R_s in ohm');
    grid(ax,'on');
end
title(layout,sprintf('Set%d: offsets around original pointwise Rs; pair lambda=%.5g Vs, iq=%.5g A',set_index,level,iq_query));
figures = fig;
if ~isempty(export_format)
    export_investigation_plot(layout,fullfile(folder,'results',sprintf('resistance_set%d',set_index)),export_format);
end
end
