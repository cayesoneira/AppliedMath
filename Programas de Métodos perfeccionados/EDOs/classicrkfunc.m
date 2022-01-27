function[x,y]=classicrkfunc(f,a,eta,h,n)
xold=a;
yold=eta;
for k=1:n
    x(k)=xold;
    y(k,:)=yold;
    
    k=f(xold,yold);
    yold=yold+(h/6)*k;
    k=f(xold+h/2,yold+(h/2)*k);
    yold=yold+(h/3)*k;
    k=f(xold+h/2,yold+(h/2)*k);
    yold=yold+(h/3)*k;
    k=f(xold+h,yold+h*k);
    yold=yold+(h/6)*k;
    
    xold=xold+h;
end
end
