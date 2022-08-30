% PRÁCTICA 6 ejercicio 7: utilizar elfin1dns.m

clear all
global h D U gamma cin L

% coeficientes de la ecuación diferencial
D=1; U=1; gamma=0.2; cin=100; L=10;
p= D;
dp= 0;
r= U;
q= gamma;
f= 0;

% extremos del intervalo
a=0; b=L;
% datos de las condiciones de contorno
datos_a=[3,-U,-U*cin];
datos_b=[2,0];
% parámetro de refinamiento de la malla
nel=10000;
l=1;
iopcoef=1;
iopblo=3;

% resolución del problema
[t,uh]=elfin1dns(p,dp,r,q,f,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo);

syms x d U Gamma c c_in L
solucion=dsolve('d*D2c-U*Dc-Gamma*c==0','-d*Dc(0)+U*c(0)==U*c_in','Dc(L)==0',x);
simplify(solucion);
solucion=subs(solucion,{d,U,Gamma,c_in,L},{1,1,0.2,100,10});
sol_exacn=double(subs(solucion,'x',t));

error=sol_exacn-uh;
norma_L2_err=norm(error,2);
disp(['La norma en L2 del error es ', num2str(norma_L2_err)])
