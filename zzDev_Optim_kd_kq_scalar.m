eval_set = load_data();
meas = eval_set.measurement;

%% Investigation settings
remove_outliers = true; % visualization only; retain the numerical factors
outlier_threshold = 6;
kernel_name = "gaussian"; % alternative: "inverse_multiquadric"
sigma = 0.4;
lambda = 1e-4;
max_centers = 1000;

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
n_centers = min(max_centers, size(Z,1));
centers = Z(randperm(size(Z,1), n_centers), :);

R2 = sum((reshape(Z, [], 1, 3) - ...
          reshape(centers, 1, [], 3)).^2, 3);

[Phi,Phi_r2,Phi_r2r2] = evaluate_rbf_kernel(R2,sigma,kernel_name);

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
    dPhi{r} = 2*(Z(:,r)-centers(:,r)').*Phi_r2;
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
        d2Phi = 4*(Z(:,r)-centers(:,r)').* ...
            (Z(:,c)-centers(:,c)').*Phi_r2r2 + 2*(r==c)*Phi_r2;
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
[Phi_ref,Phi_ref_r2,Phi_ref_r2r2] = evaluate_rbf_kernel(R2_ref,sigma,kernel_name);
dPhi_f_ref = 2*(z_f-centers(:,3)').*Phi_ref_r2;
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
% Polynomial terms are at most linear in i_f; only the RBFs contribute.
d2Phi_ff = 4*(z_f-centers(:,3)').^2.*Phi_r2r2 + 2*Phi_r2;
d2Phi_ff_ref = 4*(z_f-centers(:,3)').^2.*Phi_ref_r2r2 + 2*Phi_ref_r2;
DeltaLff = coenergy_scale/scale(3)^2 * ...
    ((d2Phi_ff-d2Phi_ff_ref)*coeff(8:end));
DeltaLff(~reference_supported) = NaN; % H; absolute Lff is still unknown

%% Lightweight diagnostics (in-sample, not independent validation)
fprintf('Kernel %s: sigma=%g, lambda=%g, centers=%d\n', ...
    kernel_name,sigma,lambda,n_centers);
flux_error = Psi-Y;
fprintf('Coenergy flux RMSE [d q] in Vs: %.6g %.6g\n', ...
    sqrt(mean(flux_error.^2,1)));
fprintf('Coenergy flux max abs error [d q] in Vs: %.6g %.6g\n', ...
    max(abs(flux_error),[],1));
fprintf('Max |Ldq-Lqd| in H: %.6g\n',max(abs(L(:,1,2)-L(:,2,1))));
fprintf('3D-hull reference (%g, %g) A: %d/%d points supported\n', ...
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

%% Visualization-only factor masks
factor_plot_mask = isfinite([kd,kq]);
for idx = 1:numel(i_f_bp)
    rows = find(i_f == i_f_bp(idx));
    for c = 1:2
        factors = [kd(rows),kq(rows)];
        finite = isfinite(factors(:,c));
        if remove_outliers && any(finite)
            center = median(factors(finite,c));
            mad_value = median(abs(factors(finite,c)-center));
            % A numerical floor handles constant/nearly constant slices.
            spread = max(1.4826*mad_value,1e-12*max(1,abs(center)));
            factor_plot_mask(rows,c) = finite & ...
                abs(factors(:,c)-center) <= outlier_threshold*spread;
        end
    end
    hidden = sum(~factor_plot_mask(rows,:),1);
    hidden_shift = sum(~all(factor_plot_mask(rows,:),2));
    fprintf('Plot mask i_f=%g A: hidden kd/kq/shifted = %d/%d/%d of %d\n', ...
        i_f_bp(idx),hidden,hidden_shift,numel(rows));
    if hidden_shift > 0.2*numel(rows)
        warning('Visualization mask hides over 20%% of the %g A slice.',i_f_bp(idx));
    end
end
shift_plot_mask = all(factor_plot_mask,2) & isfinite(s) & isfinite(g);

%% Slice geometry and field-current colors, reused by the physical plots
field_colors = turbo(numel(i_f_bp));
field_labels = compose("i_f = %.2f A",i_f_bp);
field_tags = compose("field_level_%d",1:numel(i_f_bp));
slice_rows = cell(numel(i_f_bp),1);
slice_tri = cell(numel(i_f_bp),1);
slice_support = cell(numel(i_f_bp),1);
on_slice_boundary = true(size(i_d));
for idx = 1:numel(i_f_bp)
    rows = find(i_f == i_f_bp(idx));
    % Repeated requested coordinates have the same fitted values. Keep one
    % vertex for meshing; the fit and measured scatter data are untouched.
    [xy,first,vertex_index] = unique(X(rows,1:2),'rows');
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
    if ~isempty(tri)
        slice_support{idx} = triangulation(tri,xy_scaled);
        boundary_vertices = true(size(xy,1),1);
        boundary_vertices(unique(tri)) = false;
        boundary_edges = freeBoundary(slice_support{idx});
        boundary_vertices(boundary_edges(:)) = true;
        on_slice_boundary(rows) = boundary_vertices(vertex_index);
    end
    fprintf('i_f=%g A: %d vertices, %d duplicates, %d triangles, %d rejected\n', ...
        i_f_bp(idx),size(xy,1),numel(rows)-numel(first),size(tri,1),removed);
    if isempty(tri)
        fprintf('  No reliable 2D surface: showing points for this slice.\n');
    end
end

%% Reference support in each measured field-current slice
% A 3D convex hull can bridge missing regions and different slice footprints.
% Require the reference to lie in the retained physical slice triangulation.
reference_supported_hull = reference_supported;
slice_reference_supported = false(numel(i_f_bp),1);
reference_xy = ([i_d_ref,i_q_ref]-origin(1:2))./scale(1:2);
for idx = 1:numel(i_f_bp)
    if ~isempty(slice_support{idx})
        slice_reference_supported(idx) = ...
            ~isnan(pointLocation(slice_support{idx},reference_xy));
    end
    rows = i_f == i_f_bp(idx);
    reference_supported(rows) = reference_supported_hull(rows) & ...
        slice_reference_supported(idx);
end
W_coenergy_rel(~reference_supported) = NaN;
psi_f_rel(~reference_supported) = NaN;
DeltaLff(~reference_supported) = NaN;
fprintf('Reference support: hull %d, slice-local %d/%d; newly omitted %d\n', ...
    sum(reference_supported_hull),sum(reference_supported),numel(i_f), ...
    sum(reference_supported_hull & ~reference_supported));
fprintf('Relative field flux / DeltaLff omitted points: %d / %d\n', ...
    sum(~isfinite(psi_f_rel)),sum(~isfinite(DeltaLff)));

%% Curvature diagnosis, always using every numerical value
% Same d/q energy weight: eigenvalues of this symmetric stator L block are
% meaningful. Its coenergy Hessian is 1.5 times this block. Unknown C''(i_f)
% cannot repair a negative stator principal block; full 3D eigenvalues are
% not physically identifiable without that field-only curvature.
cross_L = (Ldq+Lqd)/2;
eigen_gap = sqrt((Ldd-Lqq).^2+4*cross_L.^2);
L_stator_min = (Ldd+Lqq-eigen_gap)/2;
L_stator_max = (Ldd+Lqq+eigen_gap)/2;
negative = L_stator_min < 0;
negative_stator_points = table(i_d(negative),i_q(negative),i_f(negative), ...
    L_stator_min(negative),L_stator_max(negative),on_slice_boundary(negative), ...
    'VariableNames',{'i_d_A','i_q_A','i_f_A','min_eigen_H','max_eigen_H','on_boundary'});
negative_stator_points = sortrows(negative_stator_points,'min_eigen_H');
fprintf('Negative stator minimum eigenvalues: %d/%d; boundary/isolated: %d\n', ...
    sum(negative),numel(i_d),sum(negative & on_slice_boundary));
disp(negative_stator_points(1:min(12,height(negative_stator_points)),:));
% The full list remains in negative_stator_points; print the worst locations.
diagnostic_values = [Ldd,Lqq,Ldf,Lqf,Lfd,Lfq,DeltaLff,psi_f_rel];
diagnostic_names = ["Ldd","Lqq","Ldf","Lqf","Lfd","Lfq", ...
                    "DeltaLff (relative)","psi_f_rel (Vs)"];
diagnostic_units = [repmat("H",1,7),"Vs"];
for c = 1:size(diagnostic_values,2)
    [lowest,row] = min(diagnostic_values(:,c),[],'omitnan');
    if isnan(lowest)
        fprintf('%s: no supported finite values\n',diagnostic_names(c));
        continue
    end
    fprintf('%s: %d negative; min %.9g %s at (%.6g, %.6g, %.6g) A\n', ...
        diagnostic_names(c),sum(diagnostic_values(:,c)<0),lowest, ...
        diagnostic_units(c),X(row,:));
end
fprintf('Requested-minus-measured current RMS [d q f] in A: %.6g %.6g %.6g\n', ...
    sqrt(mean((X-[meas.id_meas(:),meas.iq_meas(:),meas.ie_meas(:)]).^2,1)));
% Negative mutual terms can reflect orientation/cross-saturation. Negative
% psi_f_rel and DeltaLff are differences and do not establish negative
% absolute field flux or Lff. Negative stator eigenvalues instead question a
% stable reversible incremental model. Smoothness/reciprocity alone do not
% constrain RBF curvature, nor prove that a measured machine is unstable.

%% Direct finite differences at supported interior operating points
% Require a neighborhood in this slice AND both adjacent field slices.
% This is an interpolation-support check, not a claim of field measurements
% between the recorded levels. Use two step sizes to expose FD sensitivity.
fd_indices = zeros(0,1);
probe_step = 2e-4*scale;
for idx = unique(round(linspace(2,numel(i_f_bp)-1,5)))
    rows = slice_rows{idx};
    rows = rows(~on_slice_boundary(rows) & reference_supported(rows));
    [~,order] = sort(sum((Z(rows,1:2)-[0.45,0.55]).^2,2));
    for candidate = rows(order)'
        probe_xy = (X(candidate,1:2) + ...
            [0 0;probe_step(1) 0;-probe_step(1) 0;0 probe_step(2);0 -probe_step(2)] ...
            - origin(1:2))./scale(1:2);
        supported = all(slice_reference_supported(idx-1:idx+1));
        for neighbor = idx-1:idx+1
            supported = supported && ~isempty(slice_support{neighbor}) && ...
                all(~isnan(pointLocation(slice_support{neighbor},probe_xy)));
        end
        if supported
            fd_indices(end+1,1) = candidate;
            break
        end
    end
end
fprintf('Finite-difference probe currents [i_d i_q i_f] in A:\n');
disp(X(fd_indices,:));
if ~isempty(fd_indices)
    relative_eval = @(query) evaluate_relative_coenergy(query,origin,scale, ...
        centers,sigma,coenergy_scale,coeff,[i_d_ref,i_q_ref],kernel_name);
    fd_expected = [Lfd(fd_indices),Lfq(fd_indices),DeltaLff(fd_indices)];
    for step_fraction = [1e-4,5e-5]
        fd_step = step_fraction*scale;
        query = X(fd_indices,:);
        [W0,~] = relative_eval(query);
        fd_field = zeros(numel(fd_indices),3);
        fd_coenergy = zeros(numel(fd_indices),3);
        for c = 1:3
            offset = zeros(1,3); offset(c) = fd_step(c);
            [Wp,fp] = relative_eval(query+offset);
            [Wm,fm] = relative_eval(query-offset);
            fd_field(:,c) = (fp-fm)/(2*fd_step(c));
            if c == 3
                fd_coenergy(:,c) = (Wp-2*W0+Wm)/fd_step(3)^2;
                fprintf('FD dW_rel/di_f max error in Vs, step %g: %.6g\n', ...
                    step_fraction,max(abs((Wp-Wm)/(2*fd_step(3))-psi_f_rel(fd_indices))));
            else
                field_offset = [0,0,fd_step(3)];
                Wpp = relative_eval(query+offset+field_offset);
                Wpm = relative_eval(query+offset-field_offset);
                Wmp = relative_eval(query-offset+field_offset);
                Wmm = relative_eval(query-offset-field_offset);
                fd_coenergy(:,c) = (Wpp-Wpm-Wmp+Wmm)/(4*fd_step(c)*fd_step(3));
            end
        end
        fprintf('FD dpsi_f_rel/d[d q f] max errors in H, step %g: %.6g %.6g %.6g\n', ...
            step_fraction,max(abs(fd_field-fd_expected),[],1));
        fprintf('FD scalar-W Hessian [fd fq ff] max errors in H: %.6g %.6g %.6g\n', ...
            max(abs(fd_coenergy-fd_expected),[],1));
    end
    fprintf('Weighted reciprocity max errors [fd fq] in H: %.6g %.6g\n', ...
        max(abs(Lfd-dq_weight*Ldf)),max(abs(Lfq-dq_weight*Lqf)));
else
    warning('No supported interior finite-difference probes were found.');
end

%% Original connectivity in shifted coordinates, only where defensible
shift_tri = cell(numel(i_f_bp),1);
for idx = 1:numel(i_f_bp)
    rows = slice_rows{idx};
    tri = slice_tri{idx};
    tri = tri(all(reshape(shift_plot_mask(rows(tri)),size(tri)),2),:);
    physical = (X(rows,1:2)-origin(1:2))./scale(1:2);
    shifted = ([s(rows),g(rows)]-origin(1:2))./scale(1:2);
    keep = false(size(tri,1),1);
    folded = false(size(tri,1),1);
    for j = 1:size(tri,1)
        vertices = tri(j,:);
        source_edges = (physical(vertices(2:3),:)-physical(vertices(1),:))';
        shifted_edges = (shifted(vertices(2:3),:)-shifted(vertices(1),:))';
        jacobian = shifted_edges/source_edges;
        stretch = svd(jacobian);
        folded(j) = det(jacobian) <= 0;
        keep(j) = ~folded(j) && min(stretch) >= 0.1 && ...
            max(stretch) <= 10 && max(stretch)/min(stretch) <= 10;
    end
    % Local checks cannot establish global injectivity. If orientation flips
    % occur anywhere in the retained slice, show points for the whole slice.
    if any(folded)
        keep(:) = false;
    end
    shift_tri{idx} = tri(keep,:);
    fprintf('Shifted i_f=%g A: %d/%d triangles retained, %d orientation flips\n', ...
        i_f_bp(idx),sum(keep),size(tri,1),sum(folded));
end

%% One figure with native tabs
fig = figure('Name','Coenergy investigation','NumberTitle','off', ...
    'WindowStyle','normal','Units','normalized','Position',[0.02 0.05 0.96 0.85]);
tabs = uitabgroup(fig);
tabFlux = uitab(tabs,'Title','Flux Linkages');
tabCoenergy = uitab(tabs,'Title','Coenergy');
tabL = uitab(tabs,'Title','Inductances');
tabK = uitab(tabs,'Title','Factors');
tabSG = uitab(tabs,'Title','s/g Flux');
tabLSG = uitab(tabs,'Title','s/g Inductances');

%% Measured stator fluxes and reconstructed field flux
t = tiledlayout(tabFlux,1,3,'TileSpacing','compact','Padding','compact');
title(t,'d/q: fitted lines and measured points | field: reconstruction only');
handles = gobjects(numel(i_f_bp),1);
psi_names = ["psi_d","psi_q"];
flux_names = [psi_names,"psi_f,rel (reconstructed)"];
flux_values = [Psi,psi_f_rel];
for c = 1:3
    ax = nexttile(t,c);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        rows = slice_rows{idx};
        h = plot_field_slice(ax,X(rows,1:2),flux_values(rows,c),slice_tri{idx}, ...
            field_colors(idx,:),field_tags(idx),field_labels(idx));
        % Keep every measurement, including repeats omitted from the mesh.
        if c <= 2
            mask = i_f == i_f_bp(idx);
            scatter3(ax,i_d(mask),i_q(mask),Y(mask,c),7,field_colors(idx,:), ...
                'filled','MarkerFaceAlpha',0.45,'MarkerEdgeAlpha',0.45, ...
                'Tag',field_tags(idx),'HandleVisibility','off');
        end
        if c == 1, handles(idx) = h; end
    end
    xlabel(ax,'i_d in A'); ylabel(ax,'i_q in A');
    zlabel(ax,flux_names(c)+" in Vs"); title(ax,flux_names(c));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Reference-relative magnetic coenergy, never absolute physical energy
t = tiledlayout(tabCoenergy,1,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(t);
hold(ax,'on');
for idx = 1:numel(i_f_bp)
    rows = slice_rows{idx};
    handles(idx) = plot_field_slice(ax,X(rows,1:2),W_coenergy_rel(rows), ...
        slice_tri{idx},field_colors(idx,:),field_tags(idx),field_labels(idx));
end
xlabel(ax,'i_d in A'); ylabel(ax,'i_q in A'); zlabel(ax,'\Delta W''_m in J');
title(ax,sprintf('Relative coenergy: reference (%g, %g) A',i_d_ref,i_q_ref));
view(ax,3); grid(ax,'on');
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

%% Fluxes in s,g: conservative wireframes and retained measured points
t = tiledlayout(tabSG,1,2,'TileSpacing','compact','Padding','compact');
title(t,'Fit: mesh where valid, otherwise points | measured: faint points');
for c = 1:2
    ax = nexttile(t,c);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        rows = slice_rows{idx};
        values = Psi(rows,c);
        values(~shift_plot_mask(rows)) = NaN;
        h = plot_field_slice(ax,[s(rows),g(rows)],values,shift_tri{idx}, ...
            field_colors(idx,:),field_tags(idx),field_labels(idx));
        mask = i_f == i_f_bp(idx) & shift_plot_mask;
        scatter3(ax,s(mask),g(mask),Y(mask,c),7,field_colors(idx,:), ...
            'filled','MarkerFaceAlpha',0.45,'MarkerEdgeAlpha',0.45, ...
            'Tag',field_tags(idx),'HandleVisibility','off');
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
        rows = slice_rows{idx};
        values = L_display(rows,panel);
        values(~shift_plot_mask(rows)) = NaN;
        h = plot_field_slice(ax,[s(rows),g(rows)],values,shift_tri{idx}, ...
            field_colors(idx,:),field_tags(idx),field_labels(idx));
        if panel == 1, handles(idx) = h; end
    end
    xlabel(ax,'s in A'); ylabel(ax,'g in A');
    zlabel(ax,L_names(panel)+" in H"); title(ax,L_names(panel));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Local factors: mask only the plotted copy, never alter kd/kq
t = tiledlayout(tabK,1,2,'TileSpacing','compact','Padding','compact');
k_values = [kd,kq];
k_names = ["k_d","k_q"];
for c = 1:2
    ax = nexttile(t,c);
    hold(ax,'on');
    for idx = 1:numel(i_f_bp)
        rows = slice_rows{idx};
        values = k_values(rows,c);
        values(~factor_plot_mask(rows,c)) = NaN;
        h = plot_field_slice(ax,X(rows,1:2),values,slice_tri{idx}, ...
            field_colors(idx,:),field_tags(idx),field_labels(idx));
        if c == 1, handles(idx) = h; end
    end
    xlabel(ax,'i_d in A'); ylabel(ax,'i_q in A');
    zlabel(ax,k_names(c)+" in A/A"); title(ax,"Local "+k_names(c));
    view(ax,3); grid(ax,'on');
end
lgd = legend(handles,field_labels,'AutoUpdate','off');
lgd.Layout.Tile = 'east';
lgd.ItemHitFcn = @toggle_field_level;

%% Small local helpers
function [phi,dphi,d2phi] = evaluate_rbf_kernel(r2,sigma,kernel_name)
% Radial derivatives are with respect to squared normalized distance r2.
switch kernel_name
    case "gaussian"
        phi = exp(-r2/sigma^2);
        dphi = -phi/sigma^2;
        d2phi = phi/sigma^4;
    case "inverse_multiquadric"
        u = 1+r2/sigma^2;
        phi = u.^(-1/2);
        dphi = -u.^(-3/2)/(2*sigma^2);
        d2phi = 3*u.^(-5/2)/(4*sigma^4);
    otherwise
        error('Use gaussian or inverse_multiquadric for kernel_name.');
end
end

function [W_rel,psi_f] = evaluate_relative_coenergy(query,origin,scale,centers, ...
    sigma,coenergy_scale,coeff,reference,kernel_name)
% Diagnostic evaluation only; callers check the query/reference support.
z = (query-origin)./scale;
zr = [repmat((reference-origin(1:2))./scale(1:2),size(z,1),1),z(:,3)];
P = [z(:,1),z(:,2),z(:,1).^2/2,z(:,1).*z(:,2),z(:,2).^2/2, ...
     z(:,1).*z(:,3),z(:,2).*z(:,3)];
Pr = [zr(:,1),zr(:,2),zr(:,1).^2/2,zr(:,1).*zr(:,2),zr(:,2).^2/2, ...
      zr(:,1).*zr(:,3),zr(:,2).*zr(:,3)];
r2 = sum((reshape(z,[],1,3)-reshape(centers,1,[],3)).^2,3);
r2_ref = sum((reshape(zr,[],1,3)-reshape(centers,1,[],3)).^2,3);
[phi,phi_r2] = evaluate_rbf_kernel(r2,sigma,kernel_name);
[phi_ref,phi_ref_r2] = evaluate_rbf_kernel(r2_ref,sigma,kernel_name);
W_rel = coenergy_scale*((P-Pr)*coeff(1:7)+(phi-phi_ref)*coeff(8:end));
psi_f = coenergy_scale/scale(3)*((z(:,1)-zr(:,1))*coeff(6) + ...
    (z(:,2)-zr(:,2))*coeff(7) + ...
    (2*(z(:,3)-centers(:,3)').*(phi_r2-phi_ref_r2))*coeff(8:end));
end

function h = plot_field_slice(ax,xy,height,tri,color,tag,label)
% Remove faces touching invalid values without retriangulating across holes.
valid = isfinite(height) & all(isfinite(xy),2);
xy(~valid,:) = NaN; % excluded points must not stretch the plotted axes
height(~valid) = NaN;
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
