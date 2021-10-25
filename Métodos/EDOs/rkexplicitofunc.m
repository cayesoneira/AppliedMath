function[x,y]=rkexplicitofunc(f,a,eta,h,n)
xold=a;
yold=eta;

c=[0,1/2,1/2,1];
A=[0,0,0,0;
    1/2,0,0,0;
    0,1/2,0,0;
    0,0,1,0];
b=[1/6,2/6,2/6,1/6];

for m=1:n
    x(m)=xold;
    y(m,:)=yold;
    k(1,:)=f(xold+c(1)*h,yold);
    for j=2:length(b)
        yi=yold;
        for i=1:j-1
            yi=yi+A(j,i)*h*k(i,:);
        end
        k(j,:)=f(xold+c(j)*h,yi);
    end
    for j=1:length(b)
        yold=yold+h*b(j)*k(j,:);
    end
    xold=xold+h;
end
end
