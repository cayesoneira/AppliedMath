hold off
clear all
clc
clf
format long

%Problema de Cauchy
a=0;
B=1;
eta=[1,1,1,1];
f=@(x,z) [2*x*z(1)*z(4),10*x*z(1)^5*z(4),2*x*z(4),-2*x*(z(3)-1)];

%Pasos
h=0.001;
n=ceil((B-a)/h);
n=n+1;
x=a:h:B;

%Solución exacta (si se da)
sol=@(x) [exp(sin(x.^2)),exp(5*sin(x.^2)),sin(x.^2)+1,cos(x.^2)];

%Funciones de los métodos particulares
tic
% [x,y]=eulerexpfunc(f,a,eta,h,n);
% [x,y]=eulerimpfunc(f,a,eta,h,n);
% [x,y]=trapeciofunc(f,a,eta,h,n);
% [x,y]=classicrkfunc(f,a,eta,h,n);
% [x,y]=rkexplicitofunc(f,a,eta,h,n);
[x,y]=ERKencajado(f,a,B,eta,h);
% [x,y]=dirkfunc(f,a,eta,h,n);
% [x,y]=gausslegendrerk2func(f,a,eta,h,n);
toc

%Error
x=x.';
soluc=sol(x);
err=abs(y-soluc);
%Da los errores máximos en ciertos nodos, no en el mismo
errormax=max(err)
errorfinal=err(end,:);

fid=1;
h=h
escribe_cabecera(fid,x(end),y(end,:),err(end,:))
escribe_paso(fid,length(x),x(end),y(end,:),err(end,:))

%Representación gráfica
if length(eta)==2
    subplot(1,3,1),plot(x,soluc(:,[1]),x,y(:,[1])),title('Solución unidimensional exacta frente a numérica')
    subplot(1,3,2),plot3(x,y(:,[1]),y(:,[2])),title('Soluciones numéricas')
    subplot(1,3,3),plot(y(:,[1]),y(:,[2])),title('Órbitas numéricas')
end
