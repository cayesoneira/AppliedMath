function[x,y]=ERKencajado(f,a,B,eta,h)
xold=a;
yold=eta;
hcur=h;

%Tolerancia que queremos conseguir del error local
eps=1.d-10;
%Orden del método principal (no del encajado)
p=2;
%k máximo y k mínimo como factores del paso y un factor
kmax=2;
kmin=0.1;
fac=0.8;

c=[0,1/2,3/4,1];
A=[0,0,0,0;
    1/2,0,0,0;
    0,3/4,0,0;
    2/9,1/3,4/9,0];
b=[2/9,1/3,4/9,0];
bt=[7/24,1/4,1/3,1/8];
m=0;
error=Inf;
kopt=0;
while xold<B
    m=m+1;
    y(m,:)=yold;
    x(m)=xold;
    error=Inf;
    while error>eps
        clear k yi
        k(1,:)=f(xold+c(1)*hcur,yold);
        for j=2:length(b)
            yi=yold;
            for i=1:j-1
                yi=yi+A(j,i)*hcur*k(i,:);
            end
            k(j,:)=f(xold+c(j)*hcur,yi);
        end
        
        errvec=zeros(length(b),1);
        for j=1:length(b)
            errvec=errvec+hcur*(b(j)-bt(j))*k(j,:);
        end
        error=norm(errvec,2);
        kopt=nthroot(eps/error,p+1);
        knew=min(kmax,max(kmin,fac*kopt));
        hnew=knew*hcur;
        hcur=hnew;
    end
    for j=1:length(b)
        yold=yold+hcur*b(j)*k(j,:);
    end
    xold=xold+hcur;
    
end
end
