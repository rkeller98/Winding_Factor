function g = fit_global_fourparameter(x,y,groups,triangles,two)
% Four coefficients, searched as [kdq kqd kd_eff kq_eff] to expose equivalence.
% Reuse the existing grid builder on pretransformed currents: model 6 with
% zero coefficients is an identity lookup in those coordinates, not a new
% common model family. The measured flux components are never transformed.
started=tic;
unit=max(abs(two.p(1)),1); center=[0 0 two.p];
scale=[.2 .2 unit/10 unit/30];
g.bounds=[-.35 -.35 two.bounds(1,:);.35 .35 two.bounds(2,:)];
g.population=two.population;
reference=four_parameter_geometry(x,groups,triangles,[0 two.p(1) 0 two.p(2)]);
% Existing measurement folds remain; do not create new folds or deepen old
% ones. Healthy triangles must retain at least 5% of their original area.
minimum_ratio=min(reference.ratio,.05);
g.seed_objective=mean(sum(two.error(g.population,:).^2,2));
% Include inward effective-factor starts: the seed can sit on active geometry
% constraints, where an all-positive initial simplex cannot move both ways.
g.start_coordinates=[center;center+[.04 -.04 -.002*unit 0];center+[-.04 .04 -.005*unit 0]];
g.solutions=[]; g.scores=[]; g.exitflags=[]; g.evaluations=[];
for j=1:size(g.start_coordinates,1)
 initial=1+(g.start_coordinates(j,:)-center)./scale;
 [v,cost,flag,out]=fminsearch(@(v) objective(center+(v-1).*scale),initial, ...
  optimset('Display','off','MaxIter',150,'MaxFunEvals',260,'TolX',2e-4,'TolFun',1e-12));
 z=center+(v-1).*scale;
 original=objective(g.start_coordinates(j,:));
 if original<cost, z=g.start_coordinates(j,:); cost=original; end
 g.solutions(j,:)=z; g.scores(j,1)=cost;
 g.exitflags(j,1)=flag; g.evaluations(j,1)=out.funcCount;
end
% Controlled comparison: optimize ONLY the shear at the unchanged two-factor
% compensation. This also supplies a geometry-equivalent admissible fallback.
[~,best]=min(g.scores);
shear_starts=unique([0 0;g.solutions(best,1:2)],'rows','stable');
shear_solutions=zeros(size(shear_starts)); shear_scores=inf(size(shear_starts,1),1);
shear_flags=zeros(size(shear_scores)); shear_evaluations=shear_flags;
for j=1:size(shear_starts,1)
 [v,cost,flag,out]=fminsearch(@(v) objective([.2*(v-1) two.p]),1+shear_starts(j,:)/.2, ...
  optimset('Display','off','MaxIter',90,'MaxFunEvals',160,'TolX',2e-4,'TolFun',1e-12));
 shear_solutions(j,:)=.2*(v-1); shear_scores(j)=cost;
 shear_flags(j)=flag; shear_evaluations(j)=out.funcCount;
end
g.shear_solutions=shear_solutions; g.shear_scores=shear_scores;
g.shear_exitflags=shear_flags; g.shear_evaluations=shear_evaluations;
[shear_cost,best]=min(shear_scores);
shear_z=[shear_solutions(best,:) two.p];
if g.seed_objective<shear_cost, shear_z=center; end
g.shear_only=solution(shear_z);
% Full overlap check on finalists, with original connectivity. An apparent
% gain involving ANY new positive-area self-overlap is not accepted.
choices=[g.solutions;shear_z;center];
costs=[g.scores;objective(shear_z);objective(center)];
[~,order]=sort(costs); rejected=0;
for j=order'
 if ~isfinite(costs(j)), continue; end
 candidate=solution(choices(j,:));
 geometry=four_parameter_geometry(x,groups,triangles,candidate.p);
 extra=0;
 for k=1:numel(triangles)
  extra=extra+nnz(~ismember(geometry.overlap_pairs{k},reference.overlap_pairs{k},'rows'));
 end
 if extra>0, rejected=rejected+1; continue; end
 break
end
g.p=candidate.p; g.effective=geometry.effective; g.geometry=geometry;
g.map=candidate.map; g.shifted=candidate.shifted;
g.prediction=candidate.prediction; g.error=candidate.error;
g.valid_rec=candidate.valid_rec; g.objective=costs(j);
g.folded=geometry.folded; g.folded_count=geometry.folded_count;
g.rejected_overlap_solutions=rejected; g.reference_geometry=reference;
% Remove the shear while retaining EXACTLY the final effective compensation.
g.effective_only.map=fit_reduction_map(x,y,6,g.effective,max(abs(x),[],1),81,true);
g.effective_only.prediction=reduction_lookup(g.effective_only.map,x);
g.effective_only.valid_rec=all(isfinite([y g.effective_only.prediction]),2);
g.effective_only.error=g.effective_only.prediction-y;
g.common=two.valid_rec & g.valid_rec;
g.control_common=g.population & g.common & g.shear_only.valid_rec & g.effective_only.valid_rec;
g.seconds=toc(started);
fprintf('  Four-parameter p=[kdq kde kqd kqe]=[%.9g %.9g %.9g %.9g], effective=[%.9g %.9g].\n',g.p,g.effective);
fprintf('  det(A)=%.8g, cond(A)=%.8g; fixed N=%d; lost=%d; supported=%d/%d; %.2fs.\n', ...
 geometry.determinant,geometry.condition,nnz(g.population),nnz(g.population & ~g.valid_rec),nnz(g.valid_rec),size(x,1),g.seconds);
fprintf('  Fixed-population mean squared vector error seed/four: %.12g / %.12g Vs^2; overlap-rejected finalists=%d.\n',g.seed_objective,g.objective,rejected);
fprintf('  Geometry two/four: reversed %d/%d, collapsed %d/%d, cond(J)>10 %d/%d, adjacent overlaps %d/%d, other self-overlaps %d/%d.\n', ...
 reference.reversed,geometry.reversed,reference.collapsed,geometry.collapsed,reference.distorted,geometry.distorted, ...
 reference.adjacent,geometry.adjacent,reference.nonadjacent,geometry.nonadjacent);
for k=1:size(g.solutions,1)
 fprintf('    Start %d: shear/effective=[%.7g %.7g %.7g %.7g], cost=%.12g, flag=%d, evaluations=%d.\n', ...
  k,g.solutions(k,:),g.scores(k),g.exitflags(k),g.evaluations(k));
end
for k=1:size(shear_solutions,1)
 fprintf('    Fixed-effective shear start %d: [%.8g %.8g], cost=%.12g, flag=%d, evaluations=%d.\n', ...
  k,shear_solutions(k,:),shear_scores(k),shear_flags(k),shear_evaluations(k));
end
predictions={two.prediction,g.prediction,g.shear_only.prediction,g.effective_only.prediction};
names={'Two-factor','Four-parameter','Shear only, original effective pair','Final effective pair, no shear'};
for k=1:4
 e=1000*(predictions{k}(g.control_common,:)-y(g.control_common,:));
 fprintf('  Control common N=%d %s RMSE [d q joint]=[%.7g %.7g %.7g] mVs.\n', ...
  nnz(g.control_common),names{k},sqrt(mean(e.^2,1)),sqrt(mean(e.^2,'all')));
end

 function p=coefficients(z)
  p=[z(1),z(3)+z(1)*z(4),z(2),z(4)+z(2)*z(3)];
 end
 function r=solution(z)
  r.p=coefficients(z);
  r.shifted=[x(:,1)+r.p(1)*x(:,2)+r.p(2)*x(:,3), ...
             r.p(3)*x(:,1)+x(:,2)+r.p(4)*x(:,3)];
  transformed=[r.shifted zeros(size(x,1),1)];
  r.map=fit_reduction_map(transformed,y,6,[0 0],[1 1 1],81,true);
  r.prediction=reduction_lookup(r.map,transformed);
  r.valid_rec=all(isfinite([r.prediction y]),2); r.error=r.prediction-y;
 end
 function cost=objective(z)
  if any(z<g.bounds(1,:) | z>g.bounds(2,:)), cost=inf; return; end
  A=[1 z(1);z(2) 1];
  if det(A)<.5 || cond(A)>3, cost=inf; return; end
  ratio=1+reference.gradient*z(3:4)';
  if any(ratio<minimum_ratio-1e-10), cost=inf; return; end
  r=solution(z);
  if any(~r.valid_rec(g.population)), cost=inf; return; end
  cost=mean(sum(r.error(g.population,:).^2,2));
 end
end
