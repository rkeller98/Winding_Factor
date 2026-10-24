function [y,J,valid] = eval_flux_model(model,x)
% Two fluxes in Vs; J(:,output,input) in Vs per physical input unit.
N=size(x,1);dim=size(x,2);J=zeros(N,2,dim);
valid=all(x>=model.bounds(1,:)-1e-10 & x<=model.bounds(2,:)+1e-10,2);
switch model.kind
 case 'rbf'
  z=(x-model.origin)./model.scale;
  dist=sum((reshape(z,N,1,dim)-reshape(model.centers,1,[],dim)).^2,3);
  phi=exp(-dist/model.width^2);y=[ones(N,1),z,phi]*model.coefficients;
  for r=1:dim
   dp=-2/model.width^2*(z(:,r)-model.centers(:,r)').*phi;
   J(:,:,r)=(model.coefficients(1+r,:)+dp*model.coefficients(dim+2:end,:))/model.scale(r);
  end
 case 'bspline'
  z=(x-model.origin)./model.scale;[A,dA]=flux_tensor_basis(z,model.knots);
  y=A*model.coefficients;
  for r=1:dim,J(:,:,r)=dA{r}*model.coefficients/model.scale(r);end
 otherwise
  map=model.map;z=reduction_coordinates(x,map.model,map.p,map.current_scale);
  z=(z-map.origin)./map.scale;nd=size(z,2);lo=zeros(N,nd);u=lo;h=lo;inside=true(N,1);
  for r=1:nd
   a=map.axes{r};ix=discretize(z(:,r),a);inside=inside&isfinite(ix);
   ix(~isfinite(ix))=1;lo(:,r)=ix;h(:,r)=reshape(a(ix+1)-a(ix),[],1);
   u(:,r)=(z(:,r)-reshape(a(ix),[],1))./h(:,r);
  end
  y=zeros(N,2);G=zeros(N,2,nd);sz=cellfun(@numel,map.axes);
  for corner=0:2^nd-1
   bits=bitget(corner,1:nd);w=ones(N,nd);ind=cell(1,nd);
   for r=1:nd,w(:,r)=(1-bits(r))*(1-u(:,r))+bits(r)*u(:,r);ind{r}=lo(:,r)+bits(r);end
   ix=sub2ind(sz,ind{:});v=[map.values{1}(ix),map.values{2}(ix)];
   y=y+v.*prod(w,2);
   for r=1:nd
    other=setdiff(1:nd,r);dw=(2*bits(r)-1)*prod(w(:,other),2)./h(:,r)/map.scale(r);
    G(:,:,r)=G(:,:,r)+v.*dw;
   end
  end
  if nd==3,J=G;else
   T=zeros(N,2,3);T(:,1,1)=1;T(:,2,2)=1;e=x(:,3);p=map.p;sc=map.current_scale;
   switch map.model
    case {1,2,3},T(:,1,3)=p(1);
    case 4,T(:,1,3)=p(1)+2*p(2)*e/sc(3);
    case 5
     T(:,1,2)=2*e*p(3).*x(:,2)/sc(2)^2;
     T(:,1,3)=p(1)+2*p(2)*e/sc(3)+p(3)*(x(:,2)/sc(2)).^2;
    case 6,T(:,1,3)=p(1);T(:,2,3)=p(2);
    case 7
     D=x(:,1)/sc(1);Q=x(:,2)/sc(2);E=e/sc(3);
     T(:,1,1)=1+e.*(p(4)+2*p(5)*D)/sc(1);
     T(:,1,2)=2*e*p(3).*Q/sc(2);
     T(:,1,3)=p(1)+2*p(2)*E+p(3)*Q.^2+p(4)*D+p(5)*D.^2;
   end
   for r=1:3,J(:,:,r)=G(:,:,1).*T(:,1,r)+G(:,:,2).*T(:,2,r);end
  end
  valid=valid&inside;
end
y(~valid,:)=nan;J(~valid,:,:)=nan;valid=valid&all(isfinite(y),2)&all(isfinite(reshape(J,N,[])),2);
end
