eval_set = load_data();
meas = eval_set.measurement;

%%

id = meas.id_req(:); 
iq = meas.iq_req(:); 
ie = meas.ie_req(:); 
psid = meas.psi_d(:); 
psiq = meas.psi_q(:); 

ie_bp = unique(ie);
ie_bp = [min(ie_bp), max(ie_bp)];
% ie_bp = ie_bp(1:2:end);


X = [id, iq, ie];
Y = [psid, psiq];

% Normalize the current coordinates
origin = min(X);
scale = max(X) - origin;

Z = (X - origin) ./ scale;

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

% Affine part + RBF basis
A = [ones(size(Z,1),1), Z, Phi];

% Fit both flux components
regularization = [zeros(n_centers,4), ...
                  sqrt(lambda)*eye(n_centers)];

W = [A; regularization] \ [Y; zeros(n_centers,2)];

L = zeros(size(Z,1), 2, 3);

for r = 1:3
    dPhi = -2/sigma^2 * ...
        (Z(:,r) - centers(:,r)').*Phi;

    L(:,:,r) = (W(r+1,:) + dPhi*W(5:end,:)) ...
               / scale(r);
end

%% Psi Plot 
Psi = A*W; 

figure()
t = tiledlayout(1,2);

axD = nexttile(t,1);
axQ = nexttile(t,2);

hold(axD, 'on');
hold(axQ, 'on');

for idx = 1:numel(ie_bp)
    mask = ie == ie_bp(idx); 
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
            mask = ie == ie_bp(idx);

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

d_mask = abs(mean(kd) - kd) > 50; 
q_mask = abs(mean(kq) - kq) > 50; 

kd(d_mask) = mean(kd);
kq(q_mask) = mean(kq); 


%% Plot Inductance Ldd and Lqq 


figure()

tiledlayout(2,2)

for idx = 1:numel(ie_bp)
    mask = ie == ie_bp(idx); 
    label = sprintf("ie = %.2f A", ie_bp(idx));
    scatter3(axD, id(mask), iq(mask), Psi(mask,1), ...
        'filled', DisplayName=label);

    scatter3(axQ, id(mask), iq(mask), Psi(mask,2), ...
        'filled', DisplayName=label);
end 



%% 

s = id + kd.*ie; 
g = iq + kq.*ie; 

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
            mask = ie == ie_bp(idx);

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
    mask = ie == ie_bp(idx); 
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