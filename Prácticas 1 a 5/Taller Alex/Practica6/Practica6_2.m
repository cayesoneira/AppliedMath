% PRÁCTICA 6 ejercicio 2

clear all
global h

%---------------------- PROBLEMA 1 ----

% coeficientes de la ecuación diferencial
fp=1; fq=-1; ff=0;

% extremos del intervalo
a=0; b=pi/2;
% datos de las condiciones de contorno
datos_a=[1,3];
datos_b=[1,7];
%parámetro de refinamiento de la malla
nel=10000;
l=1;
iopcoef=1;
iopblo=3;

% resolución del problema
[x,uh]=elfin1d(fp,fq,ff,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo);


% CÁLCULO DE LA SOLUCIÓN EXACTA
ecuacion='-D2x-x=0';
condicion='x(0)=3, x(pi/2)=7';
solucion=dsolve(ecuacion,condicion);
solucion=simplify(solucion)
u=double(subs(solucion,'t',x));

error=u-uh; N=nel+1;
error_L2_1=sqrt(h*(error(1)^2/2+sum(error(2:N-1).^2)+error(N)^2/2));
disp(['Norma L2 del error= ', num2str(error_L2_1)])


%---------------------- PROBLEMA 2----
clear all
global h

% coeficientes de la ecuación diferencial
fp2=@(t) t.^2;
fq2=@(t) 1-t; 
ff2=@(t) -t.^3-7*t.^2-2*t;

% extremos del intervalo
a=0; b=1;
% datos de las condiciones de contorno
datos_a=[1,0];
datos_b=[3,1,7];
%parámetro de refinamiento de la malla
nel=10000;
l=1;
iopcoef=2;
iopblo=3;

% resolución del problema
[x,uh]=elfin1d(fp2,fq2,ff2,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo);


% CÁLCULO DE LA SOLUCIÓN EXACTA
ecuacion='-t^2*D2x-2*t*Dx+(1-t)*x=-t^3-7*t^2-2*t';
condicion='x(0)=0, Dx(1)+x(1)=7';
solucion=dsolve(ecuacion,condicion);
solucion=simplify(solucion)
u2=double(subs(solucion,'t',x));

error=u2-uh; N=nel+1;
error_L2_2=sqrt(h*(error(1)^2/2+sum(error(2:N-1).^2)+error(N)^2/2));
disp(['Norma L2 del error= ', num2str(error_L2_2)])

