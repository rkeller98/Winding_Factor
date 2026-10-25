function g = fit_consistent_twoaxis(x,y,seed_map)
% Focused method F: EVERY candidate gets its own matching reduced LUT.
% The fit population is fixed before searching. A one-cell support guard at
% the saved method-F seed removes interpolation-boundary ambiguity; candidates
% losing ANY of these points are infeasible, not assigned a 1-Vs residual.
seed=reshape(seed_map.p,1,2);
unit=max(abs(seed(1)),1);
scale=[unit unit/4];
g.bounds=[seed-[.25*unit .1*unit];seed+[.25*unit .1*unit]];
z=(reduction_coordinates(x,6,seed,max(abs(x),[],1))-seed_map.origin)./seed_map.scale;
g.population=all(isfinite([x y]),2);
F=griddedInterpolant(seed_map.axes,seed_map.values{1},'linear','none');
hd=seed_map.axes{1}(2)-seed_map.axes{1}(1);
hq=seed_map.axes{2}(2)-seed_map.axes{2}(1);
for d=-1:1
    for q=-1:1
        g.population=g.population & isfinite(F(z(:,1)+d*hd,z(:,2)+q*hq));
    end
end
if ~any(g.population)
    error('The method-F seed has no one-cell-interior supported fit population.');
end
g.seed_p=seed;
g.seed_objective=objective(seed);
% Small, deterministic multi-start refinement of ONLY method F.
[D,Q]=ndgrid([-.05 0 .05]*unit,[-.02 0 .02]*unit);
initial=seed+[D(:),Q(:)];
scores=inf(size(D(:)));
for j=1:size(initial,1), scores(j)=objective(initial(j,:)); end
[~,order]=sort(scores);
eligible=order(isfinite(scores(order)));
g.start_pairs=unique([seed;initial(eligible(1:min(2,numel(eligible))),:)],'rows','stable');
g.solutions=nan(size(g.start_pairs));
g.scores=inf(size(g.start_pairs,1),1);
g.exitflags=zeros(size(g.scores));
g.evaluations=zeros(size(g.scores));
for j=1:size(g.start_pairs,1)
    start=1+(g.start_pairs(j,:)-seed)./scale;
    [v,cost,flag,out]=fminsearch(@(v) objective(seed+(v-1).*scale),start, ...
        optimset('Display','off','MaxIter',110,'MaxFunEvals',200,'TolX',2e-4,'TolFun',1e-12));
    p=seed+(v-1).*scale;
    original=objective(g.start_pairs(j,:));
    if original<cost, p=g.start_pairs(j,:); cost=original; end
    g.solutions(j,:)=p; g.scores(j)=cost;
    g.exitflags(j)=flag; g.evaluations(j)=out.funcCount;
end
[g.objective,best]=min(g.scores);
g.p=g.solutions(best,:);
g.map=fit_reduction_map(x,y,6,g.p,max(abs(x),[],1),81,true);
g.shifted=reduction_coordinates(x,6,g.p,max(abs(x),[],1));
g.prediction=reduction_lookup(g.map,x);
g.valid_rec=all(isfinite([g.prediction y]),2);
g.error=g.prediction-y;
g.error_mVs=1000*vecnorm(g.error,2,2);
g.initialization_spread=max(g.solutions,[],1)-min(g.solutions,[],1);
g.score_spread=max(g.scores)-min(g.scores);

    function cost=objective(p)
        if any(p<g.bounds(1,:) | p>g.bounds(2,:)), cost=inf; return; end
        candidate=fit_reduction_map(x,y,6,p,max(abs(x),[],1),81,true);
        prediction=reduction_lookup(candidate,x(g.population,:));
        if ~all(isfinite(prediction),'all'), cost=inf; return; end
        error=prediction-y(g.population,:);
        cost=mean(sum(error.^2,2));
    end
end
