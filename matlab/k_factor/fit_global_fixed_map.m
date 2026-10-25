function g = fit_global_fixed_map(x,y,map)
% Legacy fixed-map ALIGNMENT diagnostic; not a consistent reduced model.
% Its original 1-Vs missing-component penalty is retained, not silently retuned.
reference=reshape(map.p,1,2);
unit=max(abs(reference(1)),1);
scale=[unit,unit/4];
g.bounds=[reference-[unit/2,unit/4];reference+[unit/2,unit/4]];
population=all(isfinite([x y]),2);
g.population=population;
g.penalty_Vs=1; % Per missing COMPONENT, as in the existing k-model objective.
F={griddedInterpolant(map.axes,map.values{1},'linear','none'), ...
   griddedInterpolant(map.axes,map.values{2},'linear','none')};

% Restricted same-reference check for the contribution of k_q. This is the
% SAME fixed-population objective with k_q=0, not another fitted flux map.
grid=linspace(g.bounds(1,1),g.bounds(2,1),25);
score=arrayfun(@(kd) objective([kd 0]),grid);
[zero_cost,order]=sort(score);
zero_pair=[grid(order(1)),0]; zero_cost=zero_cost(1);
for j=1:3
    i=order(j);
    [kd,cost]=fminbnd(@(kd) objective([kd 0]),grid(max(1,i-1)),grid(min(numel(grid),i+1)), ...
        optimset('Display','off','TolX',1e-6,'MaxFunEvals',100));
    if cost<zero_cost, zero_pair=[kd 0]; zero_cost=cost; end
end
g.zero_q.p=zero_pair;
g.zero_q.objective=zero_cost;
g.zero_q.prediction=evaluate(zero_pair);
g.zero_q.valid_rec=population & all(isfinite(g.zero_q.prediction),2);

% Normalize around ONE rather than zero: fminsearch then creates useful
% simplex steps in both parameters. Explicit bounds reject out-of-box trials.
[D,Q]=ndgrid(linspace(g.bounds(1,1),g.bounds(2,1),9), ...
    linspace(g.bounds(1,2),g.bounds(2,2),9));
coarse=[D(:),Q(:)];
coarse_score=zeros(size(D(:)));
for j=1:size(coarse,1), coarse_score(j)=objective(coarse(j,:)); end
[~,coarse_order]=sort(coarse_score);
g.start_pairs=unique([reference;reference+[.05*unit 0];reference-[.05*unit 0]; ...
    reference+[0 .02*unit];reference-[0 .02*unit];zero_pair;coarse(coarse_order(1:3),:)], ...
    'rows','stable');
g.solutions=nan(size(g.start_pairs));
g.scores=nan(size(g.start_pairs,1),1);
g.exitflags=zeros(size(g.scores));
for j=1:size(g.start_pairs,1)
    start=1+(g.start_pairs(j,:)-reference)./scale;
    [z,cost,flag]=fminsearch(@(z) objective(reference+(z-1).*scale),start, ...
        optimset('Display','off','MaxIter',200,'MaxFunEvals',400,'TolX',1e-6,'TolFun',1e-12));
    p=reference+(z-1).*scale;
    % Retain the initial point if the non-smooth search made it worse.
    initial=objective(g.start_pairs(j,:));
    if initial<cost, p=g.start_pairs(j,:); cost=initial; end
    g.solutions(j,:)=p;
    g.scores(j)=cost;
    g.exitflags(j)=flag;
end
[g.objective,best]=min(g.scores);
g.p=g.solutions(best,:);
g.prediction=evaluate(g.p);
g.shifted=x(:,1:2)+x(:,3).*g.p;
g.valid_rec=population & all(isfinite(g.prediction),2);
g.error=g.prediction-y;
g.error_mVs=1000*vecnorm(g.error,2,2);
g.reference_prediction=evaluate(reference);
g.reference_objective=objective(reference);
g.initialization_spread=max(g.solutions,[],1)-min(g.solutions,[],1);
g.score_spread=max(g.scores)-min(g.scores);
% Compare free k_q and the restricted fit on identical supported points.
g.zero_q.common=g.valid_rec & g.zero_q.valid_rec;
e=1000*(g.prediction(g.zero_q.common,:)-y(g.zero_q.common,:));
z=1000*(g.zero_q.prediction(g.zero_q.common,:)-y(g.zero_q.common,:));
g.zero_q.joint_rmse=[sqrt(mean(z.^2,'all')),sqrt(mean(e.^2,'all'))];

    function prediction = evaluate(p)
        shifted=x(:,1:2)+x(:,3).*p;
        normalized=(shifted-map.origin)./map.scale;
        prediction=[F{1}(normalized(:,1),normalized(:,2)),F{2}(normalized(:,1),normalized(:,2))];
    end

    function cost = objective(p)
        outside=max(max(g.bounds(1,:)-p,p-g.bounds(2,:)),0);
        if any(outside>0)
            cost=1000+sum((outside./scale).^2);
            return
        end
        residual=evaluate(p)-y;
        missing=~all(isfinite(residual),2);
        residual(missing,:)=g.penalty_Vs;
        % Fixed N and equal d/q weights in Vs. Unsupported points cannot
        % disappear from the denominator; their two residuals are 1 Vs each.
        cost=mean(sum(residual(population,:).^2,2));
    end
end
