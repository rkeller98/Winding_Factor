eval_set = load_data();
meas = eval_set.measurement;

%% Measurement Data

i_d = meas.id_req(:);
i_q = meas.iq_req(:);
i_f = meas.ie_req(:);

psi_d = meas.psi_d(:);
psi_q = meas.psi_q(:);

X = [i_d, i_q, i_f];
Y = [psi_d, psi_q];

%% Physical currents and energy convention

% For this EESM, i_f denotes physical field-winding current. Requested
% currents represent the winding currents here: tracking and calibration
% are assumed, not verified. load_data applies no winding-ratio conversion;
% the underlying measurement field names stay unchanged.

% Resources/GUI_VICE-MeasurementEvalKit/SourceCode/file_prepare.m and saved
% Total.OptimParam use '3 Phase', factor 3/2. The local research notes in
% Resources/Research/magnetic_coenergy also check power and torque identities.
% Amplitude-invariant dq:
% dW_coenergy = 1.5*(psi_d*di_d + psi_q*di_q) + psi_f*di_f.
% The physical field-current orientation must match i_f. If physical field
% current = a*i_f instead, physical psi_f = (dW_coenergy/di_f)/a.
% Power-invariant dq needs factor 1. No field-turns calibration is inferred.
% A reversible, quasi-static map at fixed rotor angle/temperature is assumed;
% hysteresis, losses and operating-point drift can violate this approximation.
% Future ASM field coordinates need a defined transformation; a PMSM has
% no independently controllable physical field current. Neither is added here.
dq_weight = 3/2;

i_f_bp = unique(i_f);


% Normalize the current coordinates
origin = min(X);
scale = max(X) - origin;

Z = (X - origin) ./ scale;

%% One scalar coenergy, fitted through its stator derivatives

% Select RBF centers
rng(1)
n_centers = min(1000, size(Z,1));
centers = Z(randperm(size(Z,1), n_centers), :);

% Gaussian RBF
sigma = 0.4;
lambda = 1e-4;

R2 = sum((reshape(Z, [], 1, 3) - ...
          reshape(centers, 1, [], 3)).^2, 3);

Phi = exp(-R2 / sigma^2);

% Linear/quadratic terms give constant flux offsets and a constant Hessian.
% Omit 1, z_f, z_f^2: stator data cannot determine ANY additive C(i_f).
z_d = Z(:,1); z_q = Z(:,2); z_f = Z(:,3);
P = [z_d, z_q, z_d.^2/2, z_d.*z_q, z_q.^2/2, z_d.*z_f, z_q.*z_f];
one = ones(size(Z,1),1);
zero = zeros(size(Z,1),1);
dP = {[one, zero, z_d, z_q, zero, z_f, zero], ...
      [zero, one, zero, z_d, z_q, zero, z_f]};
dPhi = cell(1,3);
for r = 1:3
    dPhi{r} = -2/sigma^2 * (Z(:,r) - centers(:,r)').*Phi;
end

% W_coenergy = coenergy_scale * [P Phi]*coeff. This current scale (A) keeps
% coefficients in Vs and the existing ridge parameter on a flux-like scale.
coenergy_scale = dq_weight * mean(scale(1:2));
Ad = coenergy_scale/(dq_weight*scale(1)) * [dP{1}, dPhi{1}];
Aq = coenergy_scale/(dq_weight*scale(2)) * [dP{2}, dPhi{2}];
regularization = [zeros(n_centers,7), sqrt(lambda)*eye(n_centers)];
coeff = [Ad; Aq; regularization] \ [psi_d; psi_q; zeros(n_centers,1)];
% Ridge chooses one finite-basis representative, not an identified C(i_f).
W_coenergy = coenergy_scale * (P*coeff(1:7) + Phi*coeff(8:end)); % J
% W_coenergy has an arbitrary field-current-only term; only differences
% at the same i_f are observable from the stator data.
Psi = [Ad*coeff, Aq*coeff];

% Each physical Hessian entry needs 1/(scale(r)*scale(c)). Divide stator
% rows by 1.5 to get dpsi_d,q/di_c. No field self-inductance is inferred.
L = zeros(size(Z,1),2,3);
quadratic_index = [3,4,6; 4,5,7];
for r = 1:2
    for c = 1:3
        d2Phi = (4/sigma^4 * (Z(:,r)-centers(:,r)').* ...
            (Z(:,c)-centers(:,c)') - 2/sigma^2*(r==c)).*Phi;
        L(:,r,c) = coenergy_scale/(dq_weight*scale(r)*scale(c)) * ...
            (coeff(quadratic_index(r,c)) + d2Phi*coeff(8:end));
    end
end
% Weighted reciprocity: Lfd=1.5*Ldf, Lfq=1.5*Lqf, while Ldq=Lqd.
Lfd = dq_weight*L(:,1,3);
Lfq = dq_weight*L(:,2,3);

%% Relative field flux: subtract the reference curve at the SAME i_f
i_d_ref = 0; % A; near the center of the measured stator-current domain
i_q_ref = 0; % A
Z_ref = [(i_d_ref-origin(1))/scale(1)*one, ...
         (i_q_ref-origin(2))/scale(2)*one, z_f];
z_d_ref = Z_ref(:,1); z_q_ref = Z_ref(:,2);
P_ref = [z_d_ref, z_q_ref, z_d_ref.^2/2, z_d_ref.*z_q_ref, ...
         z_q_ref.^2/2, z_d_ref.*z_f, z_q_ref.*z_f];
R2_ref = sum((reshape(Z_ref,[],1,3) - reshape(centers,1,[],3)).^2,3);
Phi_ref = exp(-R2_ref/sigma^2);
dPhi_f_ref = -2/sigma^2 * (z_f-centers(:,3)').*Phi_ref;
W_coenergy_rel = coenergy_scale * ((P-P_ref)*coeff(1:7) + ...
    (Phi-Phi_ref)*coeff(8:end)); % J; W_coenergy(X)-W_coenergy(X_ref)
% Use the same expression for measured points and the reference diagnostic.
field_flux_difference = @(z_d,z_q,dPhi_f) coenergy_scale/scale(3) * ...
    ((z_d-z_d_ref)*coeff(6) + (z_q-z_q_ref)*coeff(7) + ...
     (dPhi_f-dPhi_f_ref)*coeff(8:end));
psi_f_rel = field_flux_difference(z_d,z_q,dPhi{3}); % Vs, with physical i_f
% C'(i_f) cancels exactly. Absolute psi_f and Lff remain unidentifiable.
% dpsi_f = Lfd*di_d + Lfq*di_q + Lff*di_f; at fixed i_f the last term is zero.
% L stores only the two stator rows; Lfd/Lfq above follow weighted reciprocity.
% A convex hull is only a support approximation; it cannot detect data gaps.
current_hull = delaunayTriangulation(Z);
reference_supported = ~isnan(pointLocation(current_hull,Z_ref));
W_coenergy_rel(~reference_supported) = NaN;
psi_f_rel(~reference_supported) = NaN;

% Relative field self-inductance: d^2(W_coenergy-W_coenergy_ref)/di_f^2.
% Polynomial terms are at most linear in i_f; only the Gaussians contribute.
d2Phi_ff = (4/sigma^4*(z_f-centers(:,3)').^2-2/sigma^2).*Phi;
d2Phi_ff_ref = (4/sigma^4*(z_f-centers(:,3)').^2-2/sigma^2).*Phi_ref;
DeltaLff = coenergy_scale/scale(3)^2 * ...
    ((d2Phi_ff-d2Phi_ff_ref)*coeff(8:end));
DeltaLff(~reference_supported) = NaN; % H; absolute Lff is still unknown

%% Lightweight diagnostics (in-sample, not independent validation)
flux_error = Psi-Y;
fprintf('Coenergy flux RMSE [d q] in Vs: %.6g %.6g\n', ...
    sqrt(mean(flux_error.^2,1)));
fprintf('Coenergy flux max abs error [d q] in Vs: %.6g %.6g\n', ...
    max(abs(flux_error),[],1));
fprintf('Max |Ldq-Lqd| in H: %.6g\n',max(abs(L(:,1,2)-L(:,2,1))));
fprintf('Reference (%g, %g) A: %d/%d reference-curve points supported\n', ...
    i_d_ref,i_q_ref,sum(reference_supported),numel(i_f));
psi_f_ref_rel = field_flux_difference(z_d_ref,z_q_ref,dPhi_f_ref);
fprintf('Max |psi_f_rel at reference| in Vs: %.6g (by construction)\n', ...
    max(abs(psi_f_ref_rel)));
fprintf('Differential L min/max in H: %.6g %.6g; nonfinite: %d\n', ...
    min(L,[],'all'),max(L,[],'all'),sum(~isfinite(L),'all'));
fprintf('Nonpositive Ldd/Lqq: %d / %d points (not clipped)\n', ...
    sum(L(:,1,1)<=0),sum(L(:,2,2)<=0));

% Keep the original independent RBF fit ONLY as a measurement diagnostic.
% Compare residual cost and stator curl, rather than integrating this scalar
% model twice: agreement of those paths would be automatic.
A_ind = [one,Z,Phi];
reg_ind = [zeros(n_centers,4),sqrt(lambda)*eye(n_centers)];
coeff_ind = [A_ind;reg_ind] \ [Y;zeros(n_centers,2)];
error_ind = A_ind*coeff_ind-Y;
curl_ind = (coeff_ind(3,1)+dPhi{2}*coeff_ind(5:end,1))/scale(2) - ...
           (coeff_ind(2,2)+dPhi{1}*coeff_ind(5:end,2))/scale(1);
fprintf('Independent flux RMSE [d q] in Vs: %.6g %.6g\n', ...
    sqrt(mean(error_ind.^2,1)));
fprintf('Independent Ldq-Lqd RMS/max abs in H: %.6g %.6g\n', ...
    sqrt(mean(curl_ind.^2)),max(abs(curl_ind)));
% For continuous fields on a simply connected supported slice, field-flux
% path independence requires dLfd/di_q=dLfq/di_d. For independent flux fits
% this is 1.5*d/di_f(Ldq-Lqd)=0; stator curl is a stronger energy diagnostic.
% One smooth scalar W_coenergy enforces these equalities, not measurement
% consistency. Next: measured-flux loop integrals inside each field-current
% slice and residuals by slice/repetition, then independent field-flux paths.

%% Local factors and exploratory coordinates
Ldd = L(:,1,1);
Ldq = L(:,1,2);
Ldf = L(:,1,3);
Lqd = L(:,2,1);
Lqq = L(:,2,2);
Lqf = L(:,2,3);

kd = Ldf ./ Ldd;
kq = Lqf ./ Lqq;
% Local component-cancellation slopes; cross terms need not permit a common
% reduction psi_d(s,g), psi_q(s,g). Small denominators stay visible.
fprintf('kd min/max: %.6g %.6g; kq min/max: %.6g %.6g\n', ...
    min(kd),max(kd),min(kq),max(kq));
fprintf('Nonfinite kd/kq: %d / %d\n',sum(~isfinite(kd)),sum(~isfinite(kq)));
s = i_d + kd.*i_f;
g = i_q + kq.*i_f;
% These pointwise coordinates are not a global energy-coordinate transform.
% Inductances plotted against s,g remain derivatives in physical i_d,i_q,i_f.

%% Slice geometry and field-current colors, reused by the physical plots
field_colors = turbo(numel(i_f_bp));
field_labels = compose("i_f = %.2f A",i_f_bp);
field_tags = compose("field_level_%d",1:numel(i_f_bp));
slice_rows = cell(numel(i_f_bp),1);
slice_tri = cell(numel(i_f_bp),1);
for idx = 1:numel(i_f_bp)
    rows = find(i_f == i_f_bp(idx));
    % Repeated requested coordinates have the same fitted values. Keep one
    % vertex for meshing; the fit and measured scatter data are untouched.
    [xy,first] = unique(X(rows,1:2),'rows');
    slice_rows{idx} = rows(first);
    xy_scaled = (xy-origin(1:2))./scale(1:2);
    tri = zeros(0,3);
    removed = 0;
    if size(xy,1) >= 3 && rank(xy_scaled-mean(xy_scaled,1)) == 2
        mesh = delaunayTriangulation(xy_scaled);
        tri = mesh.ConnectivityList;
        edge = edges(mesh);
        lengths = sqrt(sum((xy_scaled(edge(:,1),:)-xy_scaled(edge(:,2),:)).^2,2));
        spacing = accumarray(edge(:),[lengths;lengths],[size(xy,1),1],@min);
        % Reject edges longer than 2.5 local nearest-neighbor spacings.
        % This limits bridges across gaps; small/unresolved holes can remain.
        keep = true(size(tri,1),1);
        longest = zeros(size(tri,1),1);
        for pair = [1 2;2 3;3 1]'
            vertices = tri(:,pair);
            length_edge = sqrt(sum((xy_scaled(vertices(:,1),:)- ...
                xy_scaled(vertices(:,2),:)).^2,2));
            local_spacing = reshape(spacing(vertices),size(vertices));
            keep = keep & length_edge <= 2.5*max(local_spacing,[],2);
            longest = max(longest,length_edge);
        end
        ab = xy_scaled(tri(:,2),:)-xy_scaled(tri(:,1),:);
        ac = xy_scaled(tri(:,3),:)-xy_scaled(tri(:,1),:);
        twice_area = abs(ab(:,1).*ac(:,2)-ab(:,2).*ac(:,1));
        keep = keep & twice_area > 0.05*longest.^2; % remove thin slivers
        removed = sum(~keep);
        tri = tri(keep,:);
    end
    slice_tri{idx} = tri;
    fprintf('i_f=%g A: %d vertices, %d duplicates, %d triangles, %d rejected\n', ...
        i_f_bp(idx),size(xy,1),numel(rows)-numel(first),size(tri,1),removed);
    if isempty(tri)
        fprintf('  No reliable 2D surface: showing points for this slice.\n');
    end
end
fprintf('Relative field flux / DeltaLff omitted points: %d / %d\n', ...
    sum(~isfinite(psi_f_rel)),sum(~isfinite(DeltaLff)));

%% One figure with native tabs
fig = figure('Name','Coenergy investigation','NumberTitle','off', ...
    'WindowStyle','normal','Units','normalized','Position',[0.02 0.05 0.96 0.85]);
tabs = uitabgroup(fig);
tabPsi = uitab(tabs,'Title','Stator Flux');
tabField = uitab(tabs,'Title','Field Flux');
tabL = uitab(tabs,'Title','Inductances');
tabK = uitab(tabs,'Title','Factors');
tabSG = uitab(tabs,'Title','s/g Flux');
tabLSG = uitab(tabs,'Title','s/g Inductances');

%% Relative field flux: all field-current levels in one axes
t = tiledlayout(tabField,1,1,'TileSpacing','compact','Padding','compact');
axF = nexttile(t);
hold(axF,'on');
handles = gobjects(numel(i_f_bp),1);
for idx = 1:numel(i_f_bp)
    rows = slice_rows{idx};
    handles(idx) = plot_field_slice(axF,X(rows,1:2),psi_f_rel(rows), ...
        slice_tri{idx},field_colors(idx,:),field_tags(idx),field_labels(idx));
end
xlabel(axF,'i_d in A'); ylabel(axF,'i_q in A'); zlabel(axF,'psi_f,rel in Vs');
title(axF,sprintf('Relative to (i_d, i_q) = (%g, %g) A',i_d_ref,i_q_ref));
view(axF,3); grid(axF,'on');
lgd = legend(axF,handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Measured stator fluxes and fitted wireframes
t = tiledlayout(tabPsi,1,2,'TileSpacing','compact','Padding','compact');
title(t,'Lines: coenergy fit | Points: measurements');
psi_names = ["psi_d","psi_q"];
for c = 1:2
    ax = nexttile(t,c);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        rows = slice_rows{idx};
        h = plot_field_slice(ax,X(rows,1:2),Psi(rows,c),slice_tri{idx}, ...
            field_colors(idx,:),field_tags(idx),field_labels(idx));
        % Keep every measurement, including repeats omitted from the mesh.
        mask = i_f == i_f_bp(idx);
        scatter3(ax,i_d(mask),i_q(mask),Y(mask,c),7,field_colors(idx,:), ...
            'filled','MarkerFaceAlpha',0.45,'MarkerEdgeAlpha',0.45, ...
            'Tag',field_tags(idx),'HandleVisibility','off');
        if c == 1, handles(idx) = h; end
    end
    xlabel(ax,'i_d in A'); ylabel(ax,'i_q in A');
    zlabel(ax,psi_names(c)+" in Vs"); title(ax,psi_names(c));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Differential inductance matrix: independent cameras and z ranges
% L remains N-by-2-by-3. The third display row contains derived couplings
% and reference-relative DeltaLff, never an identified absolute Lff.
L_display = [Ldd,Ldq,Ldf,Lqd,Lqq,Lqf,Lfd,Lfq,DeltaLff];
L_names = ["Ldd","Ldq","Ldf","Lqd","Lqq","Lqf","Lfd","Lfq","\Delta Lff"];
t = tiledlayout(tabL,3,3,'TileSpacing','loose','Padding','compact');
for panel = 1:9
    ax = nexttile(t,panel);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        rows = slice_rows{idx};
        h = plot_field_slice(ax,X(rows,1:2),L_display(rows,panel), ...
            slice_tri{idx},field_colors(idx,:),field_tags(idx),field_labels(idx));
        if panel == 1, handles(idx) = h; end
    end
    xlabel(ax,'i_d in A'); ylabel(ax,'i_q in A');
    zlabel(ax,L_names(panel)+" in H"); title(ax,L_names(panel));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;
% The former separate Ldd/Lqq figure duplicated the first and fifth panels.
% Do not link absolute cameras: each inductance has a different height scale.

%% Measured fluxes in the exploratory s,g coordinates
% Use scatter here: pointwise factors can fold the coordinates, so physical
% slice triangles would create misleading surfaces after this transformation.
t = tiledlayout(tabSG,1,2,'TileSpacing','compact','Padding','compact');
for c = 1:2
    ax = nexttile(t,c);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        mask = i_f == i_f_bp(idx);
        h = scatter3(ax,s(mask),g(mask),Y(mask,c),10,field_colors(idx,:), ...
            'filled','Tag',field_tags(idx),'DisplayName',field_labels(idx));
        if c == 1, handles(idx) = h; end
    end
    xlabel(ax,'s in A'); ylabel(ax,'g in A');
    zlabel(ax,psi_names(c)+" in Vs"); title(ax,psi_names(c));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Physical stator inductances plotted against s,g
t = tiledlayout(tabLSG,2,3,'TileSpacing','loose','Padding','compact');
for panel = 1:6
    ax = nexttile(t,panel);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        mask = i_f == i_f_bp(idx);
        h = scatter3(ax,s(mask),g(mask),L_display(mask,panel),10, ...
            field_colors(idx,:),'filled','Tag',field_tags(idx), ...
            'DisplayName',field_labels(idx));
        if panel == 1, handles(idx) = h; end
    end
    xlabel(ax,'s in A'); ylabel(ax,'g in A');
    zlabel(ax,L_names(panel)+" in H"); title(ax,L_names(panel));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Local factors: points preserve large ratios near small denominators
t = tiledlayout(tabK,1,2,'TileSpacing','compact','Padding','compact');
k_values = [kd,kq];
k_names = ["k_d","k_q"];
for c = 1:2
    ax = nexttile(t,c);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        mask = i_f == i_f_bp(idx);
        h = scatter3(ax,i_d(mask),i_q(mask),k_values(mask,c),10, ...
            field_colors(idx,:),'filled','Tag',field_tags(idx), ...
            'DisplayName',field_labels(idx));
        if c == 1, handles(idx) = h; end
    end
    xlabel(ax,'i_d in A'); ylabel(ax,'i_q in A');
    zlabel(ax,k_names(c)+" in A/A"); title(ax,"Local "+k_names(c));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Small plotting helpers
function h = plot_field_slice(ax,xy,height,tri,color,tag,label)
% Remove faces touching invalid values without retriangulating across holes.
valid = isfinite(height);
tri = tri(all(reshape(valid(tri),size(tri)),2),:);
if isempty(tri)
    h = scatter3(ax,xy(valid,1),xy(valid,2),height(valid),12,color, ...
        'filled','Tag',tag,'DisplayName',label);
else
    h = trimesh(tri,xy(:,1),xy(:,2),height,'Parent',ax, ...
        'FaceColor','none','EdgeColor',color,'LineWidth',0.5, ...
        'Tag',tag,'DisplayName',label);
    used = false(size(height));
    used(tri(:)) = true;
    isolated = valid & ~used;
    if any(isolated)
        scatter3(ax,xy(isolated,1),xy(isolated,2),height(isolated),12,color, ...
            'filled','Tag',tag,'HandleVisibility','off');
    end
end
end

function toggle_field_level(lgd,event)
% findall also finds hidden objects; the legend remains clickable to restore.
% Limit the shared field-current toggle to the clicked legend's tab.
tab = ancestor(lgd,'uitab');
objects = findall(tab,'Tag',event.Peer.Tag);
if strcmp(event.Peer.Visible,'on')
    set(objects,'Visible','off');
else
    set(objects,'Visible','on');
end
end
