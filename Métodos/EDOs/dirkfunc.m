function[x,y]=dirkfunc(f,a,eta,h,n)

c=[1/4,3/4];
A=[1/4,0;
    1/2,1/4];
b=[1/2,1/2];

x(1)=a;
y(1,:)=eta;
for i=1:n-1
    o=@(t) t-f(x(i)+c(1)*h,y(i,:)+A(1,1)*h*t);
    k(1,:)=fsolve(o,f(x(i),y(i,:)),optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
    for j=2:length(b)
        v=zeros(1,length(eta));
        for l=1:(j-1)
            v=v+A(j,l)*k(l,:);
        end
        o=@(t) t-f(x(i)+c(j)*h,y(i,:)+h*(v+A(j,j)*t));
        k(j,:)=fsolve(o,f(x(i),y(i,:)),optimset('Display','off','TolFun',1.e-8,'TolX',1.e-8));
    end
    vv=zeros(1:length(eta));
    for l=1:length(b)
        vv=vv+b(l)*k(l,:);
    end
    y(i+1,:)=y(i,:)+h*vv;
    x(i+1)=x(i)+h;
end
end