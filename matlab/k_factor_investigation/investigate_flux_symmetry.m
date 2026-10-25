function [study, figures] = investigate_flux_symmetry(set_index, group_indices, export_format, data, make_plots)
% Compare mirrored physical stator currents within nominal excitation groups.
% Optional data is one canonical measurement struct, e.g. with recomputed flux.
if nargin < 1 || isempty(set_index)
    set_index = 1;
end
if nargin < 2
    group_indices = [];
end
if nargin < 3
    export_format = '';
end
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(fullfile(root, 'matlab', 'common'));
if nargin < 4 || isempty(data)
    loaded = load(fullfile(root, 'results', 'reduction_data.mat'), 'data');
    data = loaded.data{set_index};
end
if nargin < 5
    make_plots = true;
end
[x, y, groups, excitation] = reduction_inputs(data);
if isempty(group_indices)
    group_indices = ceil(numel(excitation)/2);
end
scale = max(abs(x(:,1:2)), [], 1);
scale(scale == 0) = 1;
figures = gobjects(0);
study.set = set_index;
study.group_indices = group_indices;
study.summary = nan(numel(group_indices), 9);
signed_colors = [linspace(.1,1,128)' linspace(.3,1,128)' ones(128,1); ...
    ones(128,1) linspace(1,.2,128)' linspace(1,.1,128)'];
for j = 1:numel(group_indices)
    g = group_indices(j);
    rows = groups == g;
    c = x(rows,:);
    values = [y(rows,:) c(:,3)];
    xy = c(:,1:2)./scale;
    dt = delaunayTriangulation(xy);
    tri = dt.ConnectivityList;
    edges = dt.edges;
    limit = 9*mean(sum((xy(edges(:,1),:)-xy(edges(:,2),:)).^2,2));
    good = true(size(tri,1),1);
    for edge = [1 2;1 3;2 3]'
        good = good & sum((xy(tri(:,edge(1)),:)-xy(tri(:,edge(2)),:)).^2,2) < limit;
    end
    qmax = min(max(c(:,2)), -min(c(:,2)));
    if qmax <= 0
        fprintf('Set%d group %d: no positive/negative q-current overlap; skipped.\n',set_index,g);
        continue;
    end
    [D,Q] = ndgrid(linspace(min(c(:,1)),max(c(:,1)),51),linspace(qmax/41,qmax,41));
    query = [D(:) Q(:);D(:) -Q(:)];
    [membership,bary] = pointLocation(dt,query./scale);
    supported = isfinite(membership);
    supported(supported) = good(membership(supported));
    interpolated = nan(size(query,1),3);
    for component = 1:3
        a = values(:,component);
        interpolated(supported,component) = sum(bary(supported,:).*a(tri(membership(supported),:)),2);
    end
    plus = interpolated(1:numel(D),:);
    minus = interpolated(numel(D)+1:end,:);
    common = all(isfinite([plus minus]),2);
    if ~any(common)
        fprintf('Set%d group %d: no common supported mirror queries; skipped.\n',set_index,g);
        continue;
    end
    plus(~common,:) = nan;
    minus(~common,:) = nan;
    item.nominal_excitation = excitation(g);
    item.id = D;
    item.abs_iq = Q;
    item.plus = plus;
    item.minus = minus;
    item.common = reshape(common,size(D));
    item.psi_d_even = reshape((plus(:,1)+minus(:,1))/2,size(D));
    item.psi_d_odd = reshape((plus(:,1)-minus(:,1))/2,size(D));
    item.psi_q_odd = reshape((plus(:,2)-minus(:,2))/2,size(D));
    item.psi_q_even = reshape((plus(:,2)+minus(:,2))/2,size(D));
    item.d_residual = plus(:,1)-minus(:,1);
    item.q_residual = plus(:,2)+minus(:,2);
    item.ie_difference = reshape(plus(:,3)-minus(:,3),size(D));
    study.group(j) = item;
    rms_d = sqrt(mean(item.d_residual(common).^2));
    rms_q = sqrt(mean(item.q_residual(common).^2));
    delta_ie = plus(common,3)-minus(common,3);
    study.summary(j,:) = [excitation(g),nnz(common),numel(common),rms_d,rms_q, ...
        max(abs(item.d_residual(common))),max(abs(item.q_residual(common))), ...
        sqrt(mean(delta_ie.^2)),max(abs(delta_ie))];
    fprintf('Set%d, nominal %g A: %d/%d mirror pairs; RMS d-even residual %.6g mVs, q-odd residual %.6g mVs; RMS mirror delta ie %.6g A.\n', ...
        set_index,excitation(g),nnz(common),numel(common),1000*rms_d,1000*rms_q,sqrt(mean(delta_ie.^2)));
    if make_plots
        [fig,plot_tab] = investigation_tab(set_index,'Flusssymmetrie',sprintf('%g A (Soll)',excitation(g)));
        layout = tiledlayout(plot_tab,2,3,'TileSpacing','compact');
        fields = {1000*item.psi_d_even,1000*item.psi_d_odd, ...
            1000*reshape(item.d_residual,size(D)),1000*item.psi_q_odd,1000*item.psi_q_even,item.ie_difference};
        names = {'d-flux even part (mVs)','d-flux odd part (mVs)','d mirror difference (mVs)', ...
            'q-flux odd part (mVs)','q-flux even part (mVs)','measured mirror delta i_e (A)'};
        for panel = 1:6
            ax = nexttile(layout);
            surf(ax,D,Q,fields{panel},'EdgeColor','none');
            view(ax,2);
            axis(ax,'tight');
            colorbar(ax);
            xlabel(ax,'i_d in A');
            ylabel(ax,'|i_q| in A');
            title(ax,names{panel});
            if ismember(panel,[2 3 5 6])
                colormap(ax,signed_colors);
                color_limit = max(abs(fields{panel}),[],'all','omitnan');
                color_limit = max(color_limit,eps);
                clim(ax,[-color_limit color_limit]);
            else
                colormap(ax,turbo(256));
            end
        end
        title(layout,sprintf('Set%d: nominal %g A; identical i_d and |i_q|, measured i_e may differ',set_index,excitation(g)));
        figures = fig;
        if ~isempty(export_format)
            export_investigation_plot(layout,fullfile(folder,'results',sprintf('symmetry_set%d_group%d',set_index,g)),export_format);
        end
    end
end
study.summary_columns = {'ie_ref_A','pairs','queries','d_residual_rms_Vs','q_residual_rms_Vs', ...
    'd_residual_max_Vs','q_residual_max_Vs','mirror_ie_rms_A','mirror_ie_max_A'};
fprintf('Mirrors share a nominal group, not necessarily identical measured i_e; no symmetry was enforced.\n');
fprintf('Measured ie follows each original interpolated group graph; its mirror differences limit any fixed-ie parity interpretation.\n');
end
