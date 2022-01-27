function[x,y]=eulerimpfunc(f,a,eta,h,n)
xold=a;
yold=eta;
F=@(z,y,x) [z-y-h*f(x,z)];
for k=1:n
    y(k,:)=yold;
    x(k)=xold;
    
    init=yold+h*f(xold,yold);
    
    xold=xold+h;
    g=@(x) F(x,yold,xold);
    yold=fsolve(g,init,optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
end
end