function analysis = investigate_field_reductions(set_index, local, export_format)
% All groups and original 100 levels; identical LUT construction for every method.
if nargin < 3
    export_format = '';
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab','common'),fullfile(root,'matlab','k_factor'));
loaded = load(fullfile(root,'results','reduction_data.mat'),'data');
[x,y,groups,excitation] = reduction_inputs(loaded.data{set_index});
if nargin < 2 || isempty(local)
    values = investigate_local_factors(set_index);
    local = values{1};
end
[reference_k,reference] = estimate_contour_factor(x,y,groups,false,excitation);
levels = reference.levels;
ng = numel(excitation);
locations = nan(ng,numel(levels),10,3);
geometry = cell(ng,1);
fprintf('\nSet%d: field reduction, all %d excitation groups and 100 original levels.\n',set_index,ng);
for g = 1:ng
    rows = groups==g;
    item = contour_group_statistics(x(rows,:),y(rows,1),levels);
    locations(g,:,:,:) = item.locations;
    geometry{g} = item.geometry;
end
delta_id = diff(locations(:,:,:,1),1,1);
delta_ie = diff(locations(:,:,:,3),1,1);
ratios = abs(delta_id)./delta_ie;
ratios(~isfinite(ratios) | abs(ratios)<eps) = nan;
signed = -delta_id./delta_ie;
ks = reshape(mean(ratios,[1 2],'omitnan'),1,10);
means = reshape(locations(:,:,1,:),ng,numel(levels),3);
analysis.reference_mean_difference = max(abs(means-reference.means),[],'all');
fprintf('Original edge-mean agreement: k difference %.3g, mean-current difference %.3g A.\n',ks(1)-reference_k,analysis.reference_mean_difference);
ks(1) = reference_k; % The unchanged estimator defines the reference.
fixed_k = nan(ng-1,numel(levels));
fixed_fraction = nan(size(fixed_k));
for g = 1:ng-1
    for i = 1:numel(levels)
        fixed = contour_fixed_iq([geometry{g}(i) geometry{g+1}(i)]);
        k = fixed.signed_k;
        k(~isfinite(k)) = nan;
        fixed_k(g,i) = mean(k,'omitnan');
        fixed_fraction(g,i) = mean(fixed.unique_pair);
    end
end
% A constant finite shift, not substitution of a spatial coefficient field.
component_k = nan(1,2);
component_masks = [local.valid_psi_d local.valid_psi_q];
for c = 1:2
    rows = component_masks(:,c);
    a = local.a(rows,c);
    b = local.b(rows,c);
    component_k(c) = sum(a.*b)/sum(a.^2);
end
names = [item.methods {'fixed iq (global signed mean)', ...
    'Differential psi_d (global)','Differential psi_q (global)', ...
    'Differential LS (global)','Zweiachsig (global F)'}];
parameters = [ks mean(fixed_k,'all','omitnan') component_k local.global_derivative_k local.global_F(1)];
scale = max(abs(x),[],1);
for m = 1:numel(names)
    model = 1;
    p = parameters(m);
    if m == numel(names)
        model = 6;
        p = local.global_F;
    end
    map = fit_reduction_map(x,y,model,p,scale,81,true);
    prediction = reduction_lookup(map,x);
    result.name = names{m};
    result.model = model;
    result.p = p;
    result.coordinates = reduction_coordinates(x,model,p,scale);
    result.map = map;
    result.prediction = prediction;
    result.error = prediction-y;
    result.valid = all(isfinite(prediction),2);
    analysis.methods(m) = result;
end
common = true(size(groups));
for m = 1:numel(names)
    common = common & analysis.methods(m).valid;
end
analysis.set = set_index;
analysis.x = x;
analysis.flux = y;
analysis.groups = groups;
analysis.excitation = excitation;
analysis.levels = levels;
analysis.locations = locations;
analysis.signed = signed;
analysis.ratios = ratios;
analysis.fixed_k = fixed_k;
analysis.fixed_fraction = fixed_fraction;
analysis.common = common;
analysis.rmse = nan(numel(names),2);
analysis.group_rmse = nan(ng,numel(names),2);
for m = 1:numel(names)
    e = analysis.methods(m).error;
    analysis.rmse(m,:) = 1000*sqrt(mean(e(common,:).^2,1));
    for g = 1:ng
        analysis.group_rmse(g,m,:) = 1000*sqrt(mean(e(common & groups==g,:).^2,1));
    end
end
fprintf('Field comparison: %d/%d common finite measurement points, 81x81 LUTs.\n',nnz(common),numel(common));
disp(array2table(analysis.rmse,'RowNames',names,'VariableNames',{'RMSE_d_mVs','RMSE_q_mVs'}));
plot_investigation_fields(analysis,export_format);
end
