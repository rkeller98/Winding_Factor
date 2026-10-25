function [k,detail] = estimate_contour_factor(x,y,groups)
% Sorted groups, filtered edges, middle 25--75% flux range, measured Ie.
span=max(y(:,1))-min(y(:,1));
levels=linspace(min(y(:,1))+.25*span,min(y(:,1))+.75*span,100);
glist=unique(groups); means=nan(numel(glist),100,3); centroid_error=nan(numel(glist),100);
for j=1:numel(glist)
    use=groups==glist(j); a=x(use,:); b=y(use,:);
    tri=delaunay(a(:,1),a(:,2)); edge=unique(sort([tri(:,[1 2]);tri(:,[1 3]);tri(:,[2 3])],2),'rows');
    d=sum((a(edge(:,1),1:2)-a(edge(:,2),1:2)).^2,2);
    edge=edge(d<9*mean(d),:); i=edge(:,1); h=edge(:,2);
    F=scatteredInterpolant(a(:,1:2),b(:,1),'linear','none');
    for l=1:100
        t=(levels(l)-b(i,1))./(b(h,1)-b(i,1)); hit=t>=0 & t<=1;
        means(j,l,:)=mean(a(i(hit),:)+t(hit).*(a(h(hit),:)-a(i(hit),:)),1);
        centroid_error(j,l)=F(reshape(means(j,l,1:2),1,2))-levels(l);
    end
end
signed=-diff(means(:,:,1))./diff(means(:,:,3));
ratios=abs(diff(means(:,:,1)))./diff(means(:,:,3));
ratios(~isfinite(ratios) | abs(ratios)<eps)=nan; k=mean(ratios,'all','omitnan');
detail=struct('levels',levels,'means',means,'ratios',ratios,'signed',signed, ...
    'centroid_error',centroid_error,'iq_drift',diff(means(:,:,2))./diff(means(:,:,3)));
end
