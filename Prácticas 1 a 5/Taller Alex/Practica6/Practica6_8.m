% PRÁCTICA 6 ejercicio 8

clear all
global D U Gamma cin L 

% coeficientes de la ecuación diferencial
D=1; U=1; Gamma=0.2; cin=100; L=10; p=2;
fp= @(t) -D*ones(1,length(t));
fdp= @(t) zeros(1,length(t));
fr= @(t) -U*ones(1,length(t));
fq= @(t) -Gamma*ones(1,length(t));
ff= @(t) zeros(1,length(t));


% extremos del intervalo
a=0; b=L;
% datos de las condiciones de contorno
datos_a=[3,U,U*cin];
datos_b=[2,0];
% parámetro de refinamiento de la malla
nel=100;
l=1;
iopcoef=2;
iopblo=2;

p=2;
tol=1.e-4;
maxit=200;

% resolución del problema por ELEMENTOS FINITOS (Iteración Funcional)
[xIF,phiIF]=elfin1dNoLinealIF(fp,fdp,fr,fq,ff,p,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo,tol,maxit);

% resolución del problema por DIFERENCIAS FINITAS (Iteración Funcional)
num_pas=nel;
gamma=-U/D; alfa=-U*cin/D;
delta=0; beta=0;
iopc=4;

[tIF,chIF]=diffincNoLinealIF(@u3_3IF,@v3_3IF,@w3_3IF,p,a,b,alfa,beta,num_pas,iopc,gamma,delta,tol,maxit);

error_IF=norm(phiIF'-chIF,inf);
disp(['Error absoluto máximo entre Elementos finitos y Diferencias finitas (IF): ', num2str(error_IF) ])

figure(1)
plot(xIF,phiIF), hold on
plot(tIF,chIF)
title('Iteración Funcional')
leg=legend('Elementos finitos', 'Diferencias finitas');
leg.Location='best';
hold off

% resolución del problema por ELEMENTOS FINITOS (Newton)
[xN,phiN]=elfin1dNoLinealN(fp,fdp,fr,fq,ff,p,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo,tol,maxit);

% resolución del problema por DIFERENCIAS FINITAS (Newton)
num_pas=nel;
gamma=-U/D; alfa=-U*cin/D;
delta=0; beta=0;
iopc=4;

tol=1.e-4;
maxit=1000;

[tN,chN]=diffincNoLinealN(@u3_3IF,@v3_3IF,@w3_3IF,p,a,b,alfa,beta,num_pas,iopc,gamma,delta,tol,maxit);

error_N=norm(phiN'-chN,inf);
disp(['Error absoluto máximo entre Elementos finitos y Diferencias finitas (NW): ', num2str(error_N) ])


figure(2)
plot(xN,phiN), hold on
plot(tN,chN)
title('Newton')
leg=legend('Elementos finitos', 'Diferencias finitas');
leg.Location='best';
hold off
