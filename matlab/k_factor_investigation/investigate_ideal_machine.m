function [study, figures] = investigate_ideal_machine(export_format)
% Seven small controlled examples; no measured results are modified.
% Exact reduction: s=id+k0*ie, psid=Ld*s+c*iq^2, psiq=Lq*iq+2*c*s*iq.
% Coenergy W'=Ld*s^2/2+c*s*iq^2+Lq*iq^2/2 gives both stator fluxes.
% The stator Hessian is [Ld,2*c*iq;2*c*iq,Lq+2*c*s].
% Optional export_format writes only into this investigation's results folder.
if nargin < 1
    export_format = '';
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab','common'),fullfile(root,'matlab','k_factor'));
k0 = 12;
Ld = 1e-3;
Lq = 1e-2;
excitation = [10;14];
sgrid = (-350:35:350)';
qgrid = (-300:30:300)';
dense_q = unique([(-300:60:300)';(-90:10:90)']);
names = {'uniform symmetric','symmetric unequal density','same points alternate mesh', ...
    'unequal density, c=0','unequal density, stronger curvature', ...
    'asymmetric truncation','symmetric narrower support'};
curvature = [2 2 2 0 3 2 2]*1e-6;
qsecond = {qgrid,dense_q,qgrid,dense_q,dense_q,(-180:30:300)',(-180:30:180)'};
study.k0 = k0;
study.Ld = Ld;
study.Lq = Lq;
study.coenergy = 'W = Ld*s^2/2 + c*s*iq^2 + Lq*iq^2/2, s=id+k0*ie';
study.case_names = names;
fprintf('\nIdeal exactly reducible machine: k0=%.6g, Ld=%.6g H, Lq=%.6g H.\n',k0,Ld,Lq);
fprintf('Coenergy W=0.5*Ld*s^2+c*s*iq^2+0.5*Lq*iq^2; stator reciprocity holds exactly.\n');
for h = 1:numel(names)
    c = curvature(h);
    currents = zeros(0,3);
    flux = zeros(0,2);
    groups = zeros(0,1);
    meshes = {[],[]};
    for j = 1:2
        q = qgrid;
        if j == 2
            q = qsecond{h};
        end
        [s,iq] = ndgrid(sgrid,q);
        % Small even row offsets remove rectangular-grid Delaunay ties while
        % preserving exact +iq/-iq symmetry whenever the support is symmetric.
        s = s + 1e-4*iq.^2;
        id = s-k0*excitation(j);
        ie = excitation(j)*ones(numel(s),1);
        currents = [currents;id(:),iq(:),ie]; %#ok<AGROW>
        flux = [flux;Ld*s(:)+c*iq(:).^2,Lq*iq(:)+2*c*s(:).*iq(:)]; %#ok<AGROW>
        groups = [groups;j*ones(numel(s),1)]; %#ok<AGROW>
        if h == 3 && j == 2
            % Same vertices/fluxes, choose a consistent structured cell diagonal.
            % The reference estimator still uses its original Delaunay mesh.
            % Here each cell's four psi_d samples are coplanar (two q rows):
            % flipping its diagonal changes edge sampling, not the d-flux surface.
            triangles = zeros(2*(numel(sgrid)-1)*(numel(q)-1),3);
            row = 0;
            for qi = 1:numel(q)-1
                for si = 1:numel(sgrid)-1
                    a = si+(qi-1)*numel(sgrid);
                    b = a+1;
                    d = a+numel(sgrid);
                    e = d+1;
                    row = row+1;
                    triangles(row,:) = [a b d];
                    row = row+1;
                    triangles(row,:) = [b e d];
                end
            end
            meshes{j} = triangles;
        end
    end
    % Evaluate the original 100-level estimator with no changes to its code.
    [conventional_k,reference] = estimate_contour_factor(currents,flux,groups,false,excitation);
    level_index = 50;
    level = reference.levels(level_index);
    pair = contour_pair_statistics(currents,flux,groups,[1 2],level,[],meshes);
    exact_locations = nan(2,3);
    mean_eta = nan(2,1);
    for j = 1:2
        p = pair.group(j).points;
        exact = p;
        exact(:,1) = (level-c*p(:,2).^2)/Ld-k0*p(:,3);
        exact_locations(j,:) = mean(exact,1);
        eta = Ld*(p(:,1)+k0*p(:,3))+c*p(:,2).^2-level;
        mean_eta(j) = mean(eta);
    end
    delta_ie = exact_locations(2,3)-exact_locations(1,3);
    exact_k = -(exact_locations(2,1)-exact_locations(1,1))/delta_ie;
    moment_term = (c/Ld)*diff(pair.mean_iq_squared(:,1))/delta_ie;
    residual_term = -diff(mean_eta)/(Ld*delta_ie);
    % Exact identity for edge means: k=k0+(c/Ld)*Delta(mean iq^2)/Delta ie
    %                              -Delta(mean eta)/(Ld*Delta ie).
    predicted_k = k0+moment_term+residual_term;
    s = currents(:,1)+k0*currents(:,3);
    qmax = max(abs(currents(:,2)));
    min_Lqq = Lq+2*c*min(s);
    determinant_bound = Ld*min_Lqq-(2*c*qmax)^2;
    min_eigenvalue_bound = (Ld+min_Lqq-sqrt((Ld-min_Lqq)^2+4*(2*c*qmax)^2))/2;
    fixed = pair.fixed_iq.signed_k;
    fixed = fixed(isfinite(fixed));
    fixed_range = [NaN NaN NaN];
    if ~isempty(fixed)
        fixed_range = [min(fixed),median(fixed),max(fixed)];
    end
    result.name = names{h};
    result.c = c;
    result.level_index = level_index;
    result.level = level;
    result.currents = currents;
    result.flux = flux;
    result.groups = groups;
    result.statistics = pair;
    result.conventional_global_k = conventional_k;
    result.conventional_level_k = reference.signed(1,level_index);
    result.conventional_level_absolute_k = reference.ratios(1,level_index);
    result.exact_contour_locations = exact_locations;
    result.exact_contour_k = exact_k;
    result.edge_mean_analytic_residual = mean_eta;
    result.second_moment_bias_term = moment_term;
    result.interpolation_residual_term = residual_term;
    result.identity_predicted_k = predicted_k;
    result.identity_difference = pair.signed_k(1)-predicted_k;
    result.stator_hessian_determinant_lower_bound = determinant_bound;
    result.stator_hessian_min_eigenvalue_lower_bound = min_eigenvalue_bound;
    result.stator_hessian_positive = Ld>0 && determinant_bound>0;
    result.fixed_iq_min_median_max = fixed_range;
    study.cases(h) = result;
    fprintf('\n%s; c=%.3g; level=%.8g Vs; points %d/%d.\n',names{h},c,level,sum(groups==1),sum(groups==2));
    fprintf('Original full100 absolute-numerator k %.9f (error %+.9f); original selected signed/absolute %.9f/%.9f.\n', ...
        conventional_k,conventional_k-k0,result.conventional_level_k,result.conventional_level_absolute_k);
    fprintf('Selected mesh edge signed/absolute k %.9f/%.9f; length %.9f; exact-contour same-q-weights %.9f.\n', ...
        pair.signed_k(1),pair.absolute_numerator_k(1),pair.signed_k(2),exact_k);
    fprintf('Edge means [id iq ie]: [%.6f %.6f %.6f] / [%.6f %.6f %.6f]; normalized q drift %+.9g.\n', ...
        pair.locations(1,1,:),pair.locations(2,1,:),pair.iq_drift(1));
    fprintf('Mean iq^2 %.9f/%.9f A^2; centroid flux departure %.9g/%.9g Vs.\n', ...
        pair.mean_iq_squared(:,1),pair.flux_departure(:,1));
    fprintf('Bias identity: k0 %g + second-moment term %+.9f + interpolation term %+.9f = %.9f; closure %.3g.\n', ...
        k0,moment_term,residual_term,predicted_k,result.identity_difference);
    fprintf('Fixed-q sampled min/median/max %.9f %.9f %.9f; unique %d/%d; exact analytic fixed-q k0=%g.\n', ...
        fixed_range,sum(pair.fixed_iq.unique_pair),numel(pair.fixed_iq.q),k0);
    fprintf('Queried common q interval [%.6g %.6g] A; this coverage is conditional on that interval.\n', ...
        min(pair.fixed_iq.q),max(pair.fixed_iq.q));
    fprintf('Conservative stator Hessian bounds: det >= %.9g H^2; eigenvalue >= %.9g H; positive=%d.\n', ...
        determinant_bound,min_eigenvalue_bound,result.stator_hessian_positive);
    if h == 3
        fprintf('Alternate-mesh contour crossings in original Delaunay d-flux: maximum level residual %.3g Vs.\n', ...
            max(abs(pair.group(2).delaunay_crossing_residual)));
    end
end
[figures,plot_tab] = investigation_tab('Ideale Maschine','Referenzfälle','Konturen');
layouts = gobjects(2,1);
layouts(1) = tiledlayout(plot_tab,2,2,'TileSpacing','compact');
layout = layouts(1);
selected = [1 2 3 7];
colors = lines(2);
for panel = 1:numel(selected)
    result = study.cases(selected(panel));
    pair = result.statistics;
    ax = nexttile(layout);
    hold(ax,'on');
    for j = 1:2
        geometry = pair.group(j);
        q = linspace(min(geometry.currents(:,2)),max(geometry.currents(:,2)),501)';
        id = (result.level-result.c*q.^2)/Ld-k0*excitation(j);
        plot(ax,id,q,'-','Color',colors(j,:),'LineWidth',1.1, ...
            'DisplayName',sprintf('analytic, i_e=%g A',excitation(j)));
        plot(ax,geometry.points(:,1),geometry.points(:,2),'.','Color',colors(j,:), ...
            'DisplayName',sprintf('edge intersections, group%d',j));
        loc = squeeze(pair.locations(j,1,:));
        plot(ax,loc(1),loc(2),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:), ...
            'DisplayName',sprintf('edge mean, group%d',j));
        plot(ax,result.exact_contour_locations(j,1),result.exact_contour_locations(j,2),'x', ...
            'Color',colors(j,:),'LineWidth',1.5,'HandleVisibility','off');
    end
    id0 = result.level/Ld-k0*excitation;
    plot(ax,id0,[0;0],'k+-','LineWidth',1.5,'DisplayName','true fixed-q shift');
    xlabel(ax,'physical i_d (A)');
    ylabel(ax,'i_q (A)');
    axis(ax,'equal');
    grid(ax,'on');
    title(ax,{result.name,sprintf('edge k %.4f; k_0=%g; crosses: means of exact-contour samples', ...
        pair.signed_k(1),k0)});
    legend(ax,'Location','best','FontSize',7);
end
[~,plot_tab] = investigation_tab('Ideale Maschine','Referenzfälle','Methodenvergleich und Fehlerquellen');
layouts(2) = tiledlayout(plot_tab,2,2,'TileSpacing','compact');
layout = layouts(2);
summary = nan(numel(names),6);
qmeans = nan(numel(names),2);
moments = nan(numel(names),2);
terms = nan(numel(names),2);
for h = 1:numel(names)
    result = study.cases(h);
    pair = result.statistics;
    summary(h,:) = [result.conventional_global_k,pair.signed_k(1:4)',result.fixed_iq_min_median_max(2)];
    qmeans(h,:) = pair.locations(:,1,2);
    moments(h,:) = pair.mean_iq_squared(:,1);
    terms(h,:) = [result.second_moment_bias_term,result.interpolation_residual_term];
end
ax = nexttile(layout);
plot(ax,summary,'o-');
yline(ax,k0,'k--','k_0');
ylabel(ax,'k (A/A)');
legend(ax,{'original full100 abs','selected mesh edge signed','length','median','trim','fixed-q median'}, ...
    'Location','best','FontSize',8);
title(ax,'Original estimator vs alternatives; case3 selected mesh differs');
ax = nexttile(layout);
bar(ax,qmeans);
ylabel(ax,'edge mean i_q (A)');
legend(ax,{'group1','group2'},'Location','best');
title(ax,'A zero q drift does not ensure an unbiased k');
ax = nexttile(layout);
bar(ax,moments);
ylabel(ax,'edge mean i_q^2 (A^2)');
title(ax,'Second moments depend on support and intersection weights');
ax = nexttile(layout);
bar(ax,terms,'stacked');
ylabel(ax,'signed edge k minus k_0 (A/A)');
legend(ax,{'second moment term','analytic interpolation residual term'},'Location','best');
title(ax,'Exact algebraic decomposition of the observed edge bias');
all_axes = findall(plot_tab,'Type','axes');
for i = 1:numel(all_axes)
    xticks(all_axes(i),1:numel(names));
    xticklabels(all_axes(i),names);
    xtickangle(all_axes(i),30);
    grid(all_axes(i),'on');
end
drawnow;
fprintf('\nA general mesh change can alter the interpolant and edge weights.\n');
fprintf('The controlled row-cell diagonal example has coplanar d-flux cell vertices, so its d interpolant is unchanged while edge sampling changes.\n');
fprintf('Exact-contour projected means hold the intersection q weights fixed and remove only analytic interpolation residual.\n');
fprintf('No material loss, resistance error, or physical non-reducibility is present in these examples.\n');
if ~isempty(export_format)
    for i = 1:numel(layouts)
        export_investigation_plot(layouts(i),fullfile(root,'matlab','k_factor_investigation','results', ...
            sprintf('ideal_machine_%d.png',i)),export_format);
    end
end
end
