function r = calculate_twoaxis_alignment(x,y,set_id,fallback_pair)
% Fixed method-F reference, then closest point in its bilinear flux image.
% Full-rank interior minima are exact inverse roots. Other minima lie on a
% cell edge or the cell's det(J)=0 line; all three cases are considered.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
persistent previous_inputs previous_results
if isempty(previous_inputs), previous_inputs={}; previous_results={}; end
file = fullfile(root,'results','flux_map_reduction.mat');
p=[];
reference_name='Method F';
if isfile(file)
    try
        saved=load(file,'study');
        if isfield(saved,'study') && iscell(saved.study) && numel(saved.study)>=set_id && ...
                isfield(saved.study{set_id},'p') && numel(saved.study{set_id}.p)>=6
            p=saved.study{set_id}.p{6};
        end
    catch
        p=[];
    end
end
if numel(p)~=2 || any(~isfinite(p))
    if nargin<4 || numel(fallback_pair)~=2 || any(~isfinite(fallback_pair))
        error('No usable method-F pair; provide the conventional seed pair [k_A 0].');
    end
    p=fallback_pair;
    reference_name='Conventional seed';
    warning('WindingFactor:ReferenceSeed', ...
        'Set%d: saved method F unavailable; fixed reference seed [k_A 0], lightweight two-parameter fit only.',set_id);
end
if numel(previous_inputs)>=set_id && isequaln(previous_inputs{set_id},{x,y,p,reference_name})
    r=previous_results{set_id};
    return
end
% Reuse the saved CONSTANT pair; construct the current-data LUT only once.
r.map = fit_reduction_map(x,y,6,p,max(abs(x),[],1),81,true);
r.reference_name=reference_name;
map = r.map;
axes_d = map.axes{1}; axes_q = map.axes{2};
sz = size(map.values{1});
[ii,jj] = ndgrid(1:sz(1)-1,1:sz(2)-1);
nodes = [sub2ind(sz,ii(:),jj(:)),sub2ind(sz,ii(:)+1,jj(:)), ...
    sub2ind(sz,ii(:),jj(:)+1),sub2ind(sz,ii(:)+1,jj(:)+1)];
V = [map.values{1}(:),map.values{2}(:)];
good = all(isfinite(reshape(V(nodes,:),size(nodes,1),[])),2);
nodes = nodes(good,:); ii=ii(good); jj=jj(good);
if isempty(nodes)
    error('Set%d method-F reference has no supported bilinear cells.',set_id);
end
a = V(nodes(:,1),:);
b = V(nodes(:,2),:)-a;
c = V(nodes(:,3),:)-a;
d = V(nodes(:,4),:)-a-b-c;
lo = [reshape(axes_d(ii),[],1),reshape(axes_q(jj),[],1)];
h = [reshape(axes_d(ii+1)-axes_d(ii),[],1),reshape(axes_q(jj+1)-axes_q(jj),[],1)];
lower_flux = min(cat(3,V(nodes(:,1),:),V(nodes(:,2),:),V(nodes(:,3),:),V(nodes(:,4),:)),[],3);
upper_flux = max(cat(3,V(nodes(:,1),:),V(nodes(:,2),:),V(nodes(:,3),:),V(nodes(:,4),:)),[],3);
% Unique grid edges, each with an adjacent supported cell for boundary repair.
all_edges = [nodes(:,[1 2]);nodes(:,[3 4]);nodes(:,[1 3]);nodes(:,[2 4])];
[edge,first] = unique(sort(all_edges,2),'rows');
edge_cell = mod(first-1,size(nodes,1))+1;
[ed,eq] = ind2sub(sz,edge);
edge_start = [reshape(axes_d(ed(:,1)),[],1),reshape(axes_q(eq(:,1)),[],1)];
edge_end = [reshape(axes_d(ed(:,2)),[],1),reshape(axes_q(eq(:,2)),[],1)];
edge_flux = V(edge(:,1),:);
edge_delta = V(edge(:,2),:)-edge_flux;
edge_norm = sum(edge_delta.^2,2);

% A bilinear cell's singular locus is a straight line in (u,v). Its flux
% image is quadratic, so its closest point requires only a cubic's roots.
fold = struct('cell',{},'uv',{},'step',{},'f0',{},'f1',{},'f2',{},'lower',{},'upper',{});
det0 = cross2(b,c); detu = cross2(b,d); detv = cross2(d,c);
for cell_id = 1:size(nodes,1)
    det_corners = det0(cell_id)+[0,detu(cell_id),detv(cell_id),detu(cell_id)+detv(cell_id)];
    if min(det_corners)>0 || max(det_corners)<0 || (detu(cell_id)==0 && detv(cell_id)==0)
        continue
    end
    pts = zeros(0,2);
    if detv(cell_id)~=0
        pts = [pts;0,-det0(cell_id)/detv(cell_id);1,-(det0(cell_id)+detu(cell_id))/detv(cell_id)];
    end
    if detu(cell_id)~=0
        pts = [pts;-det0(cell_id)/detu(cell_id),0;-(det0(cell_id)+detv(cell_id))/detu(cell_id),1];
    end
    pts = unique(min(max(pts(all(pts>=-1e-12 & pts<=1+1e-12,2),:),0),1),'rows');
    if size(pts,1)<2, continue; end
    [~,far] = max(sum((pts-pts(1,:)).^2,2));
    uv=pts(1,:); step=pts(far,:)-uv;
    f0=a(cell_id,:)+b(cell_id,:)*uv(1)+c(cell_id,:)*uv(2)+d(cell_id,:)*prod(uv);
    f1=b(cell_id,:)*step(1)+c(cell_id,:)*step(2)+d(cell_id,:)*(uv(1)*step(2)+uv(2)*step(1));
    f2=d(cell_id,:)*prod(step);
    probe=[f0;f0+f1+f2];
    for component=1:2
        if f2(component)~=0
            t=-f1(component)/(2*f2(component));
            if t>0 && t<1, probe(end+1,:)=f0+t*f1+t^2*f2; end
        end
    end
    fold(end+1)=struct('cell',cell_id,'uv',uv,'step',step,'f0',f0,'f1',f1,'f2',f2, ...
        'lower',min(probe,[],1),'upper',max(probe,[],1));
end

n=size(x,1);
r.factors=nan(n,2); r.shifted=nan(n,2); r.prediction=nan(n,2);
r.nonunique=false(n,1); r.ill_conditioned=false(n,1); r.unsupported=false(n,1);
r.solver_completed=false(n,1);
r.excitation_threshold=.1*median(abs(x(:,3)));
r.low_excitation=abs(x(:,3))<=r.excitation_threshold;
F={griddedInterpolant(map.axes,map.values{1},'linear','none'), ...
   griddedInterpolant(map.axes,map.values{2},'linear','none')};
for point=find(~r.low_excitation)'
    target=y(point,:);
    tolerance=128*eps(max(sum(target.^2),1e-12));
    flux_tol=sqrt(tolerance);
    t=sum((target-edge_flux).*edge_delta,2)./edge_norm;
    t(edge_norm==0)=0;
    t=min(max(t,0),1);
    error=edge_flux+t.*edge_delta-target;
    losses=sum(error.^2,2);
    [best,index]=min(losses);
    coords=edge_start(index,:)+t(index)*(edge_end(index,:)-edge_start(index,:));
    winning=edge_cell(index);
    tied=find(abs(losses-best)<=tolerance);
    alternatives=edge_start(tied,:)+t(tied).*(edge_end(tied,:)-edge_start(tied,:));
    flat_tie=any(edge_norm(tied)<=tolerance);
    candidates=find(all(target>=lower_flux-flux_tol & target<=upper_flux+flux_tol,2));
    for cell_id=candidates'
        offset=target-a(cell_id,:);
        polynomial=[-cross2(b(cell_id,:),d(cell_id,:)), ...
            cross2(offset,d(cell_id,:))-cross2(b(cell_id,:),c(cell_id,:)),cross2(offset,c(cell_id,:))];
        for u=real_roots(polynomial)'
            if u<-1e-9 || u>1+1e-9, continue; end
            u=min(max(u,0),1);
            direction=c(cell_id,:)+u*d(cell_id,:);
            if sum(direction.^2)<=realmin, continue; end
            v=sum((offset-u*b(cell_id,:)).*direction)/sum(direction.^2);
            if v<-1e-9 || v>1+1e-9, continue; end
            v=min(max(v,0),1);
            prediction=a(cell_id,:)+u*b(cell_id,:)+v*c(cell_id,:)+u*v*d(cell_id,:);
            loss=sum((prediction-target).^2);
            if loss>tolerance, continue; end
            location=lo(cell_id,:)+[u v].*h(cell_id,:);
            if loss<best-tolerance
                alternatives=location; flat_tie=false;
            elseif abs(loss-best)<=tolerance
                alternatives(end+1,:)=location;
            end
            if loss<best
                best=loss; coords=location; winning=cell_id;
            end
        end
    end
    if best>tolerance
        for f=1:numel(fold)
            part=fold(f);
            distance=max(max(part.lower-target,target-part.upper),0);
            if sum(distance.^2)>best+tolerance, continue; end
            offset=part.f0-target;
            polynomial=[2*sum(part.f2.^2),3*sum(part.f1.*part.f2), ...
                sum(part.f1.^2)+2*sum(offset.*part.f2),sum(offset.*part.f1)];
            roots_in=[0;1;real_roots(polynomial)];
            for t=roots_in(roots_in>=0 & roots_in<=1)'
                prediction=part.f0+t*part.f1+t^2*part.f2;
                loss=sum((prediction-target).^2);
                location=lo(part.cell,:)+(part.uv+t*part.step).*h(part.cell,:);
                if loss<best-tolerance
                    alternatives=location; flat_tie=false;
                elseif abs(loss-best)<=tolerance
                    alternatives(end+1,:)=location;
                end
                if loss<best
                    best=loss; coords=location; winning=part.cell;
                end
            end
        end
    end
    % These LUTs do contain multiple inverse branches: leave their factors
    % undefined. Repeated solutions at the same grid knot are one inverse.
    r.solver_completed(point)=isfinite(best);
    r.nonunique(point)=flat_tie || any(vecnorm(alternatives-coords,2,2)>1e-8);
    uv=(coords-lo(winning,:))./h(winning,:);
    J=[(b(winning,:)+uv(2)*d(winning,:))'/(h(winning,1)*map.scale(1)), ...
       (c(winning,:)+uv(1)*d(winning,:))'/(h(winning,2)*map.scale(2))];
    r.ill_conditioned(point)=~all(isfinite(J),'all') || rcond(J)<1e-6;
    if r.nonunique(point) || r.ill_conditioned(point), continue; end
    physical=map.origin+coords.*map.scale;
    factors=(physical-x(point,1:2))/x(point,3);
    shifted=x(point,1:2)+factors*x(point,3);
    normalized=(shifted-map.origin)./map.scale;
    prediction=[F{1}(normalized(1),normalized(2)),F{2}(normalized(1),normalized(2))];
    if ~all(isfinite(prediction))
        middle=map.origin+(lo(winning,:)+h(winning,:)/2).*map.scale;
        step=8*eps(max([abs(x(point,1:2));abs(physical);ones(1,2)],[],1))/abs(x(point,3)) ...
            +8*eps(max(abs(factors),1));
        factors=factors+sign((middle-physical)/x(point,3)).*step;
        shifted=x(point,1:2)+factors*x(point,3);
        normalized=(shifted-map.origin)./map.scale;
        prediction=[F{1}(normalized(1),normalized(2)),F{2}(normalized(1),normalized(2))];
    end
    inside=all(normalized>=lo(winning,:)-1e-14 & normalized<=lo(winning,:)+h(winning,:)+1e-14);
    if ~inside || ~all(isfinite(prediction)) || sum((prediction-target).^2)>best+tolerance
        r.unsupported(point)=true;
        continue
    end
    r.factors(point,:)=factors; r.shifted(point,:)=shifted; r.prediction(point,:)=prediction;
end
r.valid_rec=all(isfinite([r.factors r.prediction]),2);
r.error=r.prediction-y;
r.error_mVs=1000*vecnorm(r.error,2,2);
% Feasibility is NOT identifiability. For nonzero ie every reference-map
% coordinate is attainable by a pair, even if its optimum is nonunique.
at_zero=[F{1}((x(:,1)-map.origin(1))/map.scale(1),(x(:,2)-map.origin(2))/map.scale(2)), ...
         F{2}((x(:,1)-map.origin(1))/map.scale(1),(x(:,2)-map.origin(2))/map.scale(2))];
r.geometric_supported=all(isfinite([x y]),2) & (x(:,3)~=0 | all(isfinite(at_zero),2));
r.optimization_failed=r.geometric_supported & ~r.low_excitation & (~r.solver_completed | r.unsupported);
% Preserve the previous coverage-penalized alignment as a separate diagnostic.
r.fixed_global=fit_global_fixed_map(x,y,map);
r.fixed_global.map=map;
r.global=fit_consistent_twoaxis(x,y,map);
r.historical_cv=[];
% Reuse archived method-F holdouts only if their physical row association and
% measured-target RMSE still agree. They do NOT validate the refined model.
if exist('saved','var') && isfield(saved,'study') && numel(saved.study)>=set_id
    H=saved.study{set_id};
    if all(isfield(H,{'maps','reconstruction','metrics','recon_common'})) && ...
            numel(H.maps)>=6 && size(H.reconstruction,3)>=6 && size(H.reconstruction,1)==n && ...
            isequaln(reduction_lookup(H.maps{6},x),H.reconstruction(:,:,6))
        row=H.metrics(H.metrics(:,2)==1 & H.metrics(:,3)==6,:);
        e=1000*(H.reconstruction(H.recon_common,:,6)-y(H.recon_common,:));
        measured_rmse=[sqrt(mean(e.^2,1)),sqrt(mean(e.^2,'all'))];
        if size(row,1)==1 && max(abs(measured_rmse-row(4:6)))<1e-8
            row=H.metrics(H.metrics(:,2)==2 & H.metrics(:,3)==6,:);
            if size(row,1)==1
                r.historical_cv=struct('rmse',row(4:6),'points',row(9),'coverage',row(10));
            end
        end
    end
end
previous_inputs{set_id}={x,y,p,reference_name};
previous_results{set_id}=r;
end

function z = cross2(a,b)
z=a(:,1).*b(:,2)-a(:,2).*b(:,1);
end

function values = real_roots(p)
first=find(abs(p)>64*eps(max(max(abs(p)),realmin)),1);
if isempty(first) || first==numel(p), values=zeros(0,1); return; end
z=roots(p(first:end)/max(abs(p)));
values=real(z(abs(imag(z))<=1e-9*(1+abs(real(z)))));
end
