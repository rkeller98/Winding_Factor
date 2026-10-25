function z = reduction_coordinates(x, model, p, scale)
% One shared transformation for both flux components. Reference excitation is 0 A.
if model==8, z=x; return; end
e=x(:,3); E=e/scale(3); Q=x(:,2)/scale(2); D=x(:,1)/scale(1);
switch model
    case {1,2,3}, a=e*p(1); b=zeros(size(e));
    case 4, a=e.*(p(1)+p(2)*E); b=zeros(size(e));
    case 5, a=e.*(p(1)+p(2)*E+p(3)*Q.^2); b=zeros(size(e));
    case 6, a=e*p(1); b=e*p(2);
    case 7, a=e.*(p(1)+p(2)*E+p(3)*Q.^2+p(4)*D+p(5)*D.^2); b=zeros(size(e));
end
z=[x(:,1)+a,x(:,2)+b];
end
