function[x,y]=eulerexpfunc(f,a,eta,h,n)
xold=a;
yold=eta;
for k=1:n
    x(k)=xold;
    y(k,:)=yold;
    yold=yold+h*f(xold,yold);
    xold=xold+h;
end
end