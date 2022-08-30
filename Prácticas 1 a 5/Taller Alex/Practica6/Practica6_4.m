% PRÁCTICA 6 ejercicio 4

clear all
global h 

% coefeicientes de la ecuación diferencial
al=0.15; an=0.2; A=al*an; L=10; E=21e10; densidad=7850; g=9.8;
fp=A*E; fq=0; ff=-densidad*g*A;     % signo menos en ff porque se trabaja a compresión

% extremos del intervalo
a=0; b=L;
% datos de las condiciones de contorno
datos_a=[1,0];
datos_b=[2,25];
% parámetro de refinamiento de la malla
nel=10000;
l=1;
iopcoef=1;
iopblo=1;

% resolución del problema
[x,uh]=elfin1d(fp,fq,ff,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo);

u=-densidad*g*A/(2*A*E).*x.^2+(25+densidad*g*A*L)/(A*E).*x; % solución exacta
error=u-uh; N=nel+1;
error_L2_2=sqrt(h*(error(1)^2/2+sum(error(2:N-1).^2)+error(N)^2/2))

% MÁXIMA COMPRESIÓN
d1= @(x,h) (x(3:end)-x(1:end-2))/(2*h);     % aproximación primera derivada
figure(1)
plot(x(2:end-1), A*E*d1(uh,x(2)-x(1)))
title('Esfuerzo normal')

% en virtud de la gráfica, la máxima compresión se da en el punto x=0.
