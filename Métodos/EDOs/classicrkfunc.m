function[x,y,p]=classicrkfunc(f,a,eta,h,n)
xold=a;
yold=eta;
p=n;
for o=1:n
    if yold(1)<0
        p=o;
        break
    end
    
    if yold(1)<10^-16
        p=o;
        break
    end
    
    x(o)=xold;
    y(o,:)=yold;
    
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
