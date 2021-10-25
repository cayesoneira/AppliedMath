function[x,y,p]=gausslegendrerk2func(f,a,eta,h,n)
%Diagrama de Bucher de RK Gauss Legendre:
c=[1/2-sqrt(3)/6;1/2+sqrt(3)/6];
A=[1/4,1/4-sqrt(3)/6;
    1/4+sqrt(3)/6,1/4];
b=[1/2;1/2];

xold=a;
yold=eta;

p=n;
for o=1:n
    if yold(1)<0
        p=o;
        break
    end
    
    x(o)=xold;
    y(o,:)=yold;
    F=@(z) [z(1,:)-f(xold+h*c(1),yold+h*A(1,1)*z(1,:)+h*A(1,2)*z(2,:));...
            z(2,:)-f(xold+h*c(2),yold+h*A(2,1)*z(1,:)+h*A(2,2)*z(2,:))];
    k=fsolve(F,[f(xold,yold);f(xold,yold)],optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
    yold=yold+h*b(1)*k(1,:)+h*b(2)*k(2,:);
    xold=xold+h;
end
end