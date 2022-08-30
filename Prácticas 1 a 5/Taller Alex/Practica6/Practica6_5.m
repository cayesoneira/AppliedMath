% PRACTICA 6 ejercicio 5
% (caso no simétrico)

clear all
global h

% coeficientes de la ecuación diferencial
p= @(t) ones(1,length(t));
dp= @(t) zeros(1,length(t));
r= @(t) t;
q= @(t) 1-t;
f= @(t) -t.^3+t.^2+4*t-2;

% extremos del intervalo
a=0; b=1;
% datos de las condiciones de contorno
datos_a=[1,0];
datos_b=[3,1,7];
% parámetro de refinamiento de la malla
nel=10000;
l=1;
iopcoef=2;
iopblo=3;

% resolución del problema
[t,uh]=elfin1dns(p,dp,r,q,f,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo);

u=t.^2+2*t; % solución exacta
error=u-uh;N=length(t);
error_L2_3=sqrt(h*(error(1)^2/2+sum(error(2:N-1).^2)+error(N)^2/2))
