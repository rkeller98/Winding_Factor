function g = four_parameter_geometry(x,groups,triangles,p)
% Original per-group connectivity; overlaps mean positive-area self-overlap.
A=[1 p(1);p(3) 1];
g.determinant=det(A); g.condition=cond(A);
g.effective=(A\p([2 4])')';
s=x(:,1:2)*A'+x(:,3)*p([2 4]);
g.folded=cell(numel(triangles),1); g.overlap_pairs=g.folded;
g.gradient=[]; g.ratio=[]; g.distortion=[];
g.reversed=0; g.collapsed=0; g.adjacent=0; g.nonadjacent=0;
for k=1:numel(triangles)
 rows=find(groups==k); t=triangles{k}; z=x(rows,1:2); e=x(rows,3);
 u=z(t(:,2),:)-z(t(:,1),:); v=z(t(:,3),:)-z(t(:,1),:);
 area=u(:,1).*v(:,2)-u(:,2).*v(:,1);
 de=e(t(:,2))-e(t(:,1)); df=e(t(:,3))-e(t(:,1));
 gradient=[(de.*v(:,2)-df.*u(:,2))./area, ...
           (df.*u(:,1)-de.*v(:,1))./area];
 ratio=1+gradient*g.effective';
 shifted=s(rows,:); u=shifted(t(:,2),:)-shifted(t(:,1),:);
 v=shifted(t(:,3),:)-shifted(t(:,1),:);
 actual=(u(:,1).*v(:,2)-u(:,2).*v(:,1))./area;
 g.folded{k}=actual<=1e-8;
 g.reversed=g.reversed+nnz(actual< -1e-8);
 g.collapsed=g.collapsed+nnz(abs(actual)<=1e-8);
 g.gradient=[g.gradient;gradient]; g.ratio=[g.ratio;ratio];
 % Conditioning of the actual discrete mapping, including measured ie.
 j11=1+p(2)*gradient(:,1); j12=p(1)+p(2)*gradient(:,2);
 j21=p(3)+p(4)*gradient(:,1); j22=1+p(4)*gradient(:,2);
 norm2=j11.^2+j12.^2+j21.^2+j22.^2;
 determinant=j11.*j22-j12.*j21;
 largest=(norm2+sqrt(max(norm2.^2-4*determinant.^2,0)))/2;
 g.distortion=[g.distortion;largest./max(abs(determinant),realmin)];
 pairs=overlap_pairs(shifted,t);
 g.overlap_pairs{k}=pairs;
 if ~isempty(pairs)
  shared=zeros(size(pairs,1),1);
  for c=1:3
   shared=shared+any(t(pairs(:,1),c)==t(pairs(:,2),:),2);
  end
  g.adjacent=g.adjacent+nnz(shared>=2);
  g.nonadjacent=g.nonadjacent+nnz(shared<2);
 end
end
g.folded_count=g.reversed+g.collapsed;
g.distorted=nnz(g.distortion>10);
end

function pairs = overlap_pairs(z,t)
% Bounding boxes followed by the separating-axis test for triangle interiors.
z=(z-min(z))./(max(z)-min(z));
P=cat(3,z(t(:,1),:),z(t(:,2),:),z(t(:,3),:));
lo=min(P,[],3); hi=max(P,[],3); pairs=zeros(0,2); tolerance=1e-12;
for i=1:size(t,1)-1
 j=(i+1:size(t,1))';
 j=j(all(min(hi(i,:),hi(j,:))-max(lo(i,:),lo(j,:))>tolerance,2));
 if isempty(j), continue; end
 a=z(t(i,:),:); b1=z(t(j,1),:); b2=z(t(j,2),:); b3=z(t(j,3),:);
 keep=true(size(j));
 for edge=1:6
  if edge<=3
   d=a(mod(edge,3)+1,:)-a(edge,:);
   normal=repmat([-d(2) d(1)],numel(j),1);
  else
   B={b1,b2,b3}; c=edge-3; d=B{mod(c,3)+1}-B{c};
   normal=[-d(:,2),d(:,1)];
  end
  normal=normal./max(vecnorm(normal,2,2),realmin);
  pa=a(:,1)*normal(:,1)'+a(:,2)*normal(:,2)';
  pb=[sum(b1.*normal,2),sum(b2.*normal,2),sum(b3.*normal,2)];
  keep=keep & min(max(pa,[],1)',max(pb,[],2))-max(min(pa,[],1)',min(pb,[],2))>tolerance;
 end
 pairs=[pairs;repmat(i,nnz(keep),1),j(keep)];
end
end
