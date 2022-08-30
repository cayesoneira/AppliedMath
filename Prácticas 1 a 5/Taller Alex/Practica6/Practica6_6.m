% PRÁCTICA 6 ejercicio 6
clear all
global h

% coeficientes de la ecuación diferencial
L=10; Text=20;
p=@fpdom;
q=@fq;
f=@ff;

% extremos del intervalo
a=0; b=L;
% datos de las condiciones de contorno
datos_a=[1,50];
datos_b=[3,3,3*Text];
% parámetro de refinamiento de la malla
nel=500;
iopcoef=2;
iopblo=2;

traza='*+';
color='br';
datos_dom=[3,a,3,5,b];  % subdominio, lista de extremos
N=nel-1;
figure(1)
for l=1:2
   [x,uh]=elfind1_o2_intconblokmultidom(p,q,f,datos_dom,datos_a,datos_b,N,l,iopcoef,iopblo);
   plot(x,uh,[traza(l),color(l)])
   hold on
end


