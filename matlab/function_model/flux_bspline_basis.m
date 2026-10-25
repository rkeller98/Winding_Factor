function [B,D] = flux_bspline_basis(z,t,p)
% Cox-de Boor basis and its derivative in the normalized input coordinate.
z=z(:); last=z==t(end);z(last)=t(end)-4*eps(t(end));
N=double(z>=t(1:end-1) & z<t(2:end));lower=[];
for degree=1:p
 if degree==p,lower=N;end
 n=numel(t)-degree-1;next=zeros(numel(z),n);
 for j=1:n
  a=t(j+degree)-t(j);b=t(j+degree+1)-t(j+1);
  if a>0,next(:,j)=next(:,j)+(z-t(j))/a.*N(:,j);end
  if b>0,next(:,j)=next(:,j)+(t(j+degree+1)-z)/b.*N(:,j+1);end
 end
 N=next;
end
B=N;D=zeros(size(B));
for j=1:size(B,2)
 a=t(j+p)-t(j);b=t(j+p+1)-t(j+1);
 if a>0,D(:,j)=D(:,j)+p/a*lower(:,j);end
 if b>0,D(:,j)=D(:,j)-p/b*lower(:,j+1);end
end
end
