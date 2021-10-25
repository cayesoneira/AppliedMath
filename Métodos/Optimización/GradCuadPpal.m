clear all
clc
format long

%Damos la matriz A que define la función y el vector b:
n=20;

%Apartado a
A=zeros(n);
b=zeros(n,1);

b(1)=1;
A(1,1)=2;
for i=2:n
    A(i,i)=2;
    A(i-1,i)=-1;
    A(i,i-1)=-1;
    b(i)=1;
end

%Formulación alernativa de la función A.
% A=4*eye(n)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
% b=zeros(n,1);
% b(1)=3;
% b(end)=3*(-1)^(n+1);
% b(2:2:end-1)=-2;
% b(3:2:end-1)=2;

%Apartado b
% A=hilb(n);
% b=[137/60;
%    29/20;
%    153/140;
%    743/840;
%    1879/2520];

%Apartado c
% A=[3,0,0,-2;
%    0,3,2,0;
%    0,2,3,0;
%    -2,0,0,3];
% b=[-5;
%     0;
%     5;
%     -10];

%Paso, tolerancia e iterante inicial:
tol=1.d-6;
eps=1.d-15; %Número máquina aprox.
itmax=3000;
xold=20:-1:1;
xold=xold.';
xold

%Ejecución
% tic
% [x,NormaDelGradiente,PasoDeConvengencia]=GradCuad(A,b,xold,itmax,tol)
% [x,NormaDelGradiente,PasoDeConvengencia]=GradCuadOptimizado(A,b,xold,itmax,tol)
% [x,NormaDelGradiente,PasoDeConvengencia]=GradConjCuad(A,b,xold,itmax,tol,eps)
[x,NormaDelGradiente,PasoDeConvengencia,normas]=GradConjCuadOptimizado(A,b,xold,itmax,tol,eps);
toc

NormaDelGradiente=NormaDelGradiente
normas=normas.'
x(3)
x(20)