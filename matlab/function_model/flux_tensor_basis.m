function [A,dA] = flux_tensor_basis(z,knots)
% Sparse tensor products, retaining only nonzero one-dimensional bases.
N=size(z,1);dim=size(z,2);index=ones(N,1);value=ones(N,1);deriv=zeros(N,1,dim);stride=1;
for r=1:dim
 [B,D]=flux_bspline_basis(z(:,r),knots{r},3);n=size(B,2);
 active=zeros(N,4);v=active;dv=active;
 for i=1:N
  ix=find(B(i,:)~=0 | D(i,:)~=0);active(i,1:numel(ix))=ix;
  v(i,1:numel(ix))=B(i,ix);dv(i,1:numel(ix))=D(i,ix);
 end
 oldvalue=value;oldindex=index;oldd=deriv;k=size(value,2);
 value=zeros(N,4*k);index=ones(N,4*k);deriv=zeros(N,4*k,dim);
 for h=1:4
  cols=(h-1)*k+(1:k);index(:,cols)=oldindex+max(active(:,h)-1,0)*stride;
  value(:,cols)=oldvalue.*v(:,h);
  for q=1:dim
   if q==r,deriv(:,cols,q)=oldvalue.*dv(:,h);else,deriv(:,cols,q)=oldd(:,:,q).*v(:,h);end
  end
 end
 stride=stride*n;
end
rows=repmat((1:N)',1,size(index,2));
A=sparse(rows(:),index(:),value(:),N,stride);
if nargout>1
 dA=cell(1,dim);
 for r=1:dim,v=deriv(:,:,r);dA{r}=sparse(rows(:),index(:),v(:),N,stride);end
end
end
