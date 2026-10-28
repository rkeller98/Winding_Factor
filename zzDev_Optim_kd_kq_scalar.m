eval_set = load_data();
meas = eval_set.measurement;

%% Physical currents and energy convention

% load_data maps these channels directly from PE_Id/Iq_meas_avg and
% I_exc_meas, without a winding-ratio conversion. Use measured currents for
% physical derivatives; requested excitation only identifies plotting groups.
id = meas.id_meas(:);
iq = meas.iq_meas(:);
ie = meas.ie_meas(:);
ie_group = meas.ie_req(:);
psid = meas.psi_d(:); 
psiq = meas.psi_q(:); 

% Resources/GUI_VICE-MeasurementEvalKit/SourceCode/file_prepare.m and saved
% Total.OptimParam use '3 Phase', factor 3/2. The local research notes in
% Resources/Research/magnetic_coenergy also check power and torque identities.
% Amplitude-invariant dq: dW = 1.5*(psi_d*did + psi_q*diq) + psi_r*die.
% Assume I_exc_meas is the physical rotor winding current ir, with the same
% positive orientation. Sensor calibration/rotor turns are not established
% here. If ir=a*ie, psi_r=(dW/die)/a instead. Power-invariant dq needs factor 1.
% A reversible, quasi-static map at fixed rotor angle/temperature is assumed;
% hysteresis, losses and operating-point drift can violate this approximation.
dq_weight = 3/2;

ie_bp = unique(ie_group);
ie_bp = [min(ie_bp), max(ie_bp)];
% ie_bp = ie_bp(1:2:end);


X = [id, iq, ie];
Y = [psid, psiq];

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
% Omit 1, ze, ze^2: stator data cannot determine ANY additive C(ie).
zd = Z(:,1); zq = Z(:,2); ze = Z(:,3);
P = [zd, zq, zd.^2/2, zd.*zq, zq.^2/2, zd.*ze, zq.*ze];
one = ones(size(Z,1),1);
zero = zeros(size(Z,1),1);
dP = {[one, zero, zd, zq, zero, ze, zero], ...
      [zero, one, zero, zd, zq, zero, ze]};
dPhi = cell(1,3);
for r = 1:3
    dPhi{r} = -2/sigma^2 * (Z(:,r) - centers(:,r)').*Phi;
end

% W = coenergy_scale * [P Phi]*coeff. This current scale (A) keeps
% coefficients in Vs and the existing ridge parameter on a flux-like scale.
coenergy_scale = dq_weight * mean(scale(1:2));
Ad = coenergy_scale/(dq_weight*scale(1)) * [dP{1}, dPhi{1}];
Aq = coenergy_scale/(dq_weight*scale(2)) * [dP{2}, dPhi{2}];
regularization = [zeros(n_centers,7), sqrt(lambda)*eye(n_centers)];
coeff = [Ad; Aq; regularization] \ [psid; psiq; zeros(n_centers,1)];
% Ridge chooses one finite-basis representative, not an identified C(ie).
Psi = [Ad*coeff, Aq*coeff];

% Each physical Hessian entry needs 1/(scale(r)*scale(c)). Divide stator
% rows by 1.5 to get dpsi_d,q/di_c. No rotor self-inductance is inferred.
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
% Weighted reciprocity: Lrd=1.5*Lde, Lrq=1.5*Lqe, while Ldq=Lqd.
Lrd = dq_weight*L(:,1,3);
Lrq = dq_weight*L(:,2,3);

%% Relative rotor flux: subtract the reference curve at the SAME ie
id_ref = 0; % A; near the center of the measured stator-current domain
iq_ref = 0; % A
Z_ref = [(id_ref-origin(1))/scale(1)*one, ...
         (iq_ref-origin(2))/scale(2)*one, ze];
zd_ref = Z_ref(:,1); zq_ref = Z_ref(:,2);
P_ref = [zd_ref, zq_ref, zd_ref.^2/2, zd_ref.*zq_ref, ...
         zq_ref.^2/2, zd_ref.*ze, zq_ref.*ze];
R2_ref = sum((reshape(Z_ref,[],1,3) - reshape(centers,1,[],3)).^2,3);
Phi_ref = exp(-R2_ref/sigma^2);
dPhi_e_ref = -2/sigma^2 * (ze-centers(:,3)').*Phi_ref;
W_relative = coenergy_scale * ((P-P_ref)*coeff(1:7) + ...
    (Phi-Phi_ref)*coeff(8:end)); % J; W(id,iq,ie)-W(id_ref,iq_ref,ie)
% Use the same expression for measured points and the reference diagnostic.
rotor_flux_difference = @(z_d,z_q,dPhi_e) coenergy_scale/scale(3) * ...
    ((z_d-zd_ref)*coeff(6) + (z_q-zq_ref)*coeff(7) + ...
     (dPhi_e-dPhi_e_ref)*coeff(8:end));
psi_r_rel = rotor_flux_difference(zd,zq,dPhi{3}); % Vs, assuming ie=ir
% C'(ie) cancels exactly. Absolute psi_r and Lrr remain unidentifiable.
% A convex hull is only a support approximation; it cannot detect data gaps.
current_hull = delaunayTriangulation(Z);
reference_supported = ~isnan(pointLocation(current_hull,Z_ref));
W_relative(~reference_supported) = NaN;
psi_r_rel(~reference_supported) = NaN;

%% Lightweight diagnostics (in-sample, not independent validation)
flux_error = Psi-Y;
fprintf('Coenergy flux RMSE [d q] in Vs: %.6g %.6g\n', ...
    sqrt(mean(flux_error.^2,1)));
fprintf('Coenergy flux max abs error [d q] in Vs: %.6g %.6g\n', ...
    max(abs(flux_error),[],1));
fprintf('Max |Ldq-Lqd| in H: %.6g\n',max(abs(L(:,1,2)-L(:,2,1))));
fprintf('Reference (%g, %g) A: %d/%d reference-curve points supported\n', ...
    id_ref,iq_ref,sum(reference_supported),numel(ie));
psi_r_ref_rel = rotor_flux_difference(zd_ref,zq_ref,dPhi_e_ref);
fprintf('Max |psi_r_rel at reference| in Vs: %.6g (by construction)\n', ...
    max(abs(psi_r_ref_rel)));
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
% For continuous fields on a simply connected supported slice, rotor path
% independence requires dLrd/diq=dLrq/did. For independent flux fits this is
% 1.5*d/die(Ldq-Lqd)=0; stator curl itself is a stronger energy diagnostic.
% One smooth scalar W enforces these equalities, not measurement consistency.
% Next: small measured-flux loop integrals inside each excitation slice and
% residuals by slice/repetition; only then compare independent rotor paths.

%% Relative rotor flux plot
figure()
t = tiledlayout(1,numel(ie_bp));
for idx = 1:numel(ie_bp)
    axR = nexttile(t,idx);
    axR.Toolbar.Visible = 'off';
    mask = ie_group == ie_bp(idx) & reference_supported;
    scatter3(axR,id(mask),iq(mask),psi_r_rel(mask),12,'filled');
    xlabel(axR,'i_d in A');
    ylabel(axR,'i_q in A');
    zlabel(axR,'psi_r,rel in Vs');
    title(axR,sprintf('i_e group = %.2f A',ie_bp(idx)));
    grid(axR,'on');
    view(axR,3);
end
% Points use measured ie within each requested group, not exact constant-ie
% surfaces. Unsupported reference points are omitted, not extrapolated.

%% Psi Plot 

figure()
t = tiledlayout(1,2);

axD = nexttile(t,1);
axQ = nexttile(t,2);

hold(axD, 'on');
hold(axQ, 'on');

for idx = 1:numel(ie_bp)
    mask = ie_group == ie_bp(idx);
    label = sprintf("ie = %.2f A", ie_bp(idx));
    scatter3(axD, id(mask), iq(mask), Psi(mask,1), ...
        'filled', DisplayName=label);

    scatter3(axQ, id(mask), iq(mask), Psi(mask,2), ...
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

L_name = ["Ldd", "Ldq", "Lde";
          "Lqd", "Lqq", "Lqe"];

ax = gobjects(2,3);

for r = 1:2
    for c = 1:3
        ax(r,c) = nexttile(t, (r-1)*3+c);
        hold(ax(r,c), 'on');

        for idx = 1:numel(ie_bp)
            mask = ie_group == ie_bp(idx);

            scatter3(ax(r,c), id(mask), iq(mask), L(mask,r,c), ...
                12, 'filled', ...
                DisplayName=sprintf('i_e = %.2f A', ie_bp(idx)));
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
Lde = L(:,1,3);
Lqd = L(:,2,1);
Lqq = L(:,2,2);
Lqe = L(:,2,3);

kd = Lde ./ Ldd;
kq = Lqe ./ Lqq;

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

for idx = 1:numel(ie_bp)
    mask = ie_group == ie_bp(idx);
    label = sprintf("ie = %.2f A", ie_bp(idx));
    scatter3(axDd, id(mask), iq(mask), Ldd(mask), ...
        'filled', DisplayName=label);

    scatter3(axQq, id(mask), iq(mask), Lqq(mask), ...
        'filled', DisplayName=label);
end 
xlabel(axDd,'i_d in A'); ylabel(axDd,'i_q in A'); zlabel(axDd,'Ldd in H');
xlabel(axQq,'i_d in A'); ylabel(axQq,'i_q in A'); zlabel(axQq,'Lqq in H');
grid(axDd,'on'); grid(axQq,'on');
view(axDd,3); view(axQq,3);
legend(axDd,'show'); legend(axQq,'show');

%% 

s = id + kd.*ie; 
g = iq + kq.*ie; 
% Keep these definitions as exploratory coordinates. Pointwise kd/kq do
% not define a global energy-coordinate transform. L below still contains
% derivatives with respect to physical id,iq,ie, only plotted against s,g.

figure()
tiledlayout(1,2)

nexttile
scatter3(s,g,psid)
xlabel("s"), ylabel("g")
nexttile
scatter3(s,g,psiq)
xlabel("s"), ylabel("g")

%% Inductance Plot 

figure()
t = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

L_name = ["Ldd", "Ldq", "Lde";
          "Lqd", "Lqq", "Lqe"];

ax = gobjects(2,3);

for r = 1:2
    for c = 1:3
        ax(r,c) = nexttile(t, (r-1)*3+c);
        hold(ax(r,c), 'on');

        for idx = 1:numel(ie_bp)
            mask = ie_group == ie_bp(idx);

            scatter3(ax(r,c), s(mask), g(mask), L(mask,r,c), ...
                12, 'filled', ...
                DisplayName=sprintf('i_e = %.2f A', ie_bp(idx)));
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

for idx = 1:numel(ie_bp)
    mask = ie_group == ie_bp(idx);
    L_name = sprintf("ie = %.2f A", ie_bp(idx));
    scatter3(axD, id(mask), iq(mask), kd(mask), ...
        'filled', DisplayName=L_name);

    scatter3(axQ, id(mask), iq(mask), kq(mask), ...
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
