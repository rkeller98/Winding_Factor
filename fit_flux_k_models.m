function [p,info,sc] = fit_flux_k_models(x,y,g,n)
sc=max(abs(x),[],1); k=estimate_contour_factor(x,y,g); unit=max(abs(k),1);
ug=unique(g); inner=ug(2:end-1); folds=unique(inner(round(linspace(1,numel(inner),min(3,numel(inner))))));
parts=cell(numel(folds),1);
for j=1:numel(folds)
    tr=g~=folds(j); ix=find(g==folds(j)); ix=ix(round(linspace(1,numel(ix),min(120,numel(ix)))));
    base=fit_reduction_map(x(tr,:),y(tr,:),1,k,sc,n,true);
    keep=all(isfinite(reduction_lookup(base,x(ix,:))),2); ix=ix(keep);
    parts{j}=struct('x',x(tr,:),'y',y(tr,:),'query',x(ix,:),'target',y(ix,:));
end
p=cell(8,1); p{1}=k; p{8}=[]; info=cell(8,1); info{1}=struct('inner_groups',folds);
for m=2:3
    fun=@(v) objective(v*unit,m,parts,sc,n,k);
    seeds=linspace(-.5,2.5,9); score=arrayfun(fun,seeds); [~,order]=sort(score);
    best=inf; chosen=nan; total=0;
    for j=1:2
        h=order(j); lo=seeds(max(1,h-1)); hi=seeds(min(numel(seeds),h+1));
        [v,fval,exitflag,output]=fminbnd(fun,lo,hi,optimset('TolX',.002,'Display','off'));
        total=total+output.funcCount;
        if fval<best,best=fval; chosen=v; flag=exitflag;end
    end
    p{m}=chosen*unit; info{m}=struct('objective_Vs2',best,'exitflag',flag, ...
        'evaluations',total+9,'inner_groups',folds,'search_range',[-.5 2.5]*unit);
end
count=[1 1 1 2 3 2 5];
for m=4:7
    base=zeros(1,count(m)); base(1)=p{3}/unit;
    % Nonzero search coordinates give every coefficient the same initial step.
    % A zero-valued physical coefficient must not get MATLAB's tiny zero step.
    start=ones(1,count(m));
    fun=@(v) objective((base+v-1)*unit,m,parts,sc,n,k);
    [v,fval,flag,output]=fminsearch(fun,start, ...
        optimset('MaxIter',220,'MaxFunEvals',350,'TolX',.003,'TolFun',1e-10,'Display','off'));
    p{m}=(base+v-1)*unit; info{m}=struct('objective_Vs2',fval,'exitflag',flag, ...
        'evaluations',output.funcCount,'inner_groups',folds);
end
end

function cost = objective(p,m,parts,sc,n,k)
unit=max(abs(k),1); v=p/unit;
if any(abs(v)>4),cost=1+sum(v.^2);return;end
if m==7
    [D,E]=ndgrid([-sc(1) sc(1)],[-sc(3) sc(3)]);
    jac=1+E(:).*(p(4)/sc(1)+2*p(5)*D(:)/sc(1)^2);
    if min(jac)<.25,cost=1+sum((jac-.25).^2);return;end
end
cost=0;
for j=1:numel(parts)
    d=parts{j}; map=fit_reduction_map(d.x,d.y,m,p,sc,n,true);
    pred=reduction_lookup(map,d.query); r=pred-d.target;
    if m==2,r=r(:,1);end
    r(~isfinite(r))=1; % Fixed denominator: lost support cannot improve the fit.
    cost=cost+mean(r(:).^2)/numel(parts);
end
if m>=4, cost=cost+1e-8*sum(v(2:end).^2); end
end
