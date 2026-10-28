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
i_f_bp = [min(i_f_bp), max(i_f_bp)];
% i_f_bp = i_f_bp(1:2:end);


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

%% Relative field flux plot
figure()
t = tiledlayout(1,numel(i_f_bp));
for idx = 1:numel(i_f_bp)
    axF = nexttile(t,idx);
    axF.Toolbar.Visible = 'off';
    mask = i_f == i_f_bp(idx) & reference_supported;
    scatter3(axF,i_d(mask),i_q(mask),psi_f_rel(mask),12,'filled');
    xlabel(axF,'i_d in A');
    ylabel(axF,'i_q in A');
    zlabel(axF,'psi_f,rel in Vs');
    title(axF,sprintf('i_f = %.2f A',i_f_bp(idx)));
    grid(axF,'on');
    view(axF,3);
end
% Slices use requested i_f; physical currents may differ through tracking.
% Unsupported reference points are omitted, not extrapolated.

%% Psi Plot 

figure()
t = tiledlayout(1,2);

axD = nexttile(t,1);
axQ = nexttile(t,2);

hold(axD, 'on');
hold(axQ, 'on');

for idx = 1:numel(i_f_bp)
    mask = i_f == i_f_bp(idx);
    label = sprintf("i_f = %.2f A", i_f_bp(idx));
    scatter3(axD, i_d(mask), i_q(mask), Psi(mask,1), ...
        'filled', DisplayName=label);

    scatter3(axQ, i_d(mask), i_q(mask), Psi(mask,2), ...
        'filled', DisplayName=label);
end 

xlabel(axD, 'i_d in A');
ylabel(axD, 'i_q in A');
zlabel(axD, 'psi_d in Vs');

xlabel(axQ, 'i_d in A');
ylabel(axQ, 'i_q in A');
zlabel(axQ, 'psi_q in Vs');

view(axD, 3);
view(axQ, 3);

grid(axD, 'on');
grid(axQ, 'on');

legend(axD, 'show');
legend(axQ, 'show');

%% Inductance Plot 

figure()
t = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

L_name = ["Ldd", "Ldq", "Ldf";
          "Lqd", "Lqq", "Lqf"];

ax = gobjects(2,3);

for r = 1:2
    for c = 1:3
        ax(r,c) = nexttile(t, (r-1)*3+c);
        hold(ax(r,c), 'on');

        for idx = 1:numel(i_f_bp)
            mask = i_f == i_f_bp(idx);

            scatter3(ax(r,c), i_d(mask), i_q(mask), L(mask,r,c), ...
                12, 'filled', ...
                DisplayName=sprintf('i_f = %.2f A', i_f_bp(idx)));
        end

        xlabel(ax(r,c), 'i_d in A');
        ylabel(ax(r,c), 'i_q in A');
        zlabel(ax(r,c), L_name(r,c) + " in H");
        title(ax(r,c), L_name(r,c));

        grid(ax(r,c), 'on');
        view(ax(r,c), 3);
        legend(ax(r,c), 'show');
    end
end

linkprop(ax(:), {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'});



%%

Ldd = L(:,1,1);
Ldq = L(:,1,2);
Ldf = L(:,1,3);
Lqd = L(:,2,1);
Lqq = L(:,2,2);
Lqf = L(:,2,3);

kd = Ldf ./ Ldd;
kq = Lqf ./ Lqq;

% These original ratios are local component-cancellation slopes: a common
% reduction psi_d(s,g), psi_q(s,g) need not exist when cross terms matter.
% Small/zero Ldd or Lqq can make the ratios large/nonfinite; leave visible.
fprintf('kd min/max: %.6g %.6g; kq min/max: %.6g %.6g\n', ...
    min(kd),max(kd),min(kq),max(kq));
fprintf('Nonfinite kd/kq: %d / %d\n',sum(~isfinite(kd)),sum(~isfinite(kq)));


%% Plot Inductance Ldd and Lqq 


figure()
t = tiledlayout(1,2);
axDd = nexttile(t,1);
axQq = nexttile(t,2);
hold(axDd,'on');
hold(axQq,'on');

for idx = 1:numel(i_f_bp)
    mask = i_f == i_f_bp(idx);
    label = sprintf("i_f = %.2f A", i_f_bp(idx));
    scatter3(axDd, i_d(mask), i_q(mask), Ldd(mask), ...
        'filled', DisplayName=label);

    scatter3(axQq, i_d(mask), i_q(mask), Lqq(mask), ...
        'filled', DisplayName=label);
end 
xlabel(axDd,'i_d in A'); ylabel(axDd,'i_q in A'); zlabel(axDd,'Ldd in H');
xlabel(axQq,'i_d in A'); ylabel(axQq,'i_q in A'); zlabel(axQq,'Lqq in H');
grid(axDd,'on'); grid(axQq,'on');
view(axDd,3); view(axQq,3);
legend(axDd,'show'); legend(axQq,'show');

%% 

s = i_d + kd.*i_f;
g = i_q + kq.*i_f;
% Keep these definitions as exploratory coordinates. Pointwise kd/kq do
% not define a global energy-coordinate transform. L below still contains
% derivatives with respect to physical i_d,i_q,i_f, only plotted against s,g.

figure()
tiledlayout(1,2)

nexttile
scatter3(s,g,psi_d)
xlabel("s in A"), ylabel("g in A"), zlabel("psi_d in Vs")
nexttile
scatter3(s,g,psi_q)
xlabel("s in A"), ylabel("g in A"), zlabel("psi_q in Vs")

%% Inductance Plot 

figure()
t = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

L_name = ["Ldd", "Ldq", "Ldf";
          "Lqd", "Lqq", "Lqf"];

ax = gobjects(2,3);

for r = 1:2
    for c = 1:3
        ax(r,c) = nexttile(t, (r-1)*3+c);
        hold(ax(r,c), 'on');

        for idx = 1:numel(i_f_bp)
            mask = i_f == i_f_bp(idx);

            scatter3(ax(r,c), s(mask), g(mask), L(mask,r,c), ...
                12, 'filled', ...
                DisplayName=sprintf('i_f = %.2f A', i_f_bp(idx)));
        end

        xlabel(ax(r,c), 's in A');
        ylabel(ax(r,c), 'g in A');
        zlabel(ax(r,c), L_name(r,c) + " in H");
        title(ax(r,c), L_name(r,c));

        grid(ax(r,c), 'on');
        view(ax(r,c), 3);
        legend(ax(r,c), 'show');
    end
end

linkprop(ax(:), {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'});



%% Plot kd kq 
figure()
t = tiledlayout(1,2);

axD = nexttile(t,1);
axQ = nexttile(t,2);

hold(axD, 'on');
hold(axQ, 'on');

for idx = 1:numel(i_f_bp)
    mask = i_f == i_f_bp(idx);
    L_name = sprintf("i_f = %.2f A", i_f_bp(idx));
    scatter3(axD, i_d(mask), i_q(mask), kd(mask), ...
        'filled', DisplayName=L_name);

    scatter3(axQ, i_d(mask), i_q(mask), kq(mask), ...
        'filled', DisplayName=L_name);
end 

xlabel(axD, 'i_d in A');
ylabel(axD, 'i_q in A');
zlabel(axD, 'k_d in A/A');
title(axD, 'Local k_d');

xlabel(axQ, 'i_d in A');
ylabel(axQ, 'i_q in A');
zlabel(axQ, 'k_q in A/A');
title(axQ, 'Local k_q');

view(axD, 3);
view(axQ, 3);

grid(axD, 'on');
grid(axQ, 'on');

legend(axD, 'show');
legend(axQ, 'show');
