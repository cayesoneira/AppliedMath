function [x,uh]=elfind1_o2_nel_d(p,q,f,a,b,alfa,beta,nel)
% [x,uh]=elfind1_o2(p,q,f,a,b,alfa,beta,gb,N)
% elfind1_o2 calcula la solución discreta del problema 
% de contorno lineal de segundo orden:
%****************************************************
% 				-(p u'(x))'+q u(x)=f  en (a,b)
% 				u(a)=alfa;	u(b)=beta, 
% Se suponen los coeficientes p,q,f, alfa y beta constantes
%****************************************************
% utilizando elementos finitos de Lagrange de grado k=1.
% nel número de elementos de la malla.
%****************************************************
% Construcción de la malla 1-dimensional.
%****************************************************
global h
k=1;%orden del espacio de elementos finitos.
Nh=nel+1
nver=nel+1; % numero de vértices: a_0,a_1,...,a_N+1
h=(b-a)/nel; % tamaño de la discretizacion
M=k*nel+1-2;
nodos=M+2;%número total de nodos
i=1:nodos; x=a+(i-1)*h/k;
%i : puntero de nodos; 
%x : vector de coordenadas de los nodos
%***********************************************
% Inicialización de la matriz A
%***********************************************
A=spalloc(nodos,nodos,3*nodos-2);
% Inicialización de c
%************************************************
c=zeros(nodos,1);
%*************************************************
% Bucle en elementos
%*************************************************
for kk=1:nel
	%*************************************************
	% Calculo matriz elemental
	%*************************************************
   Kk=p/h*[1, -1;-1, 1];Mk=(q*h)/6*[2 1;1 2];Ak=Mk+Kk;
	%*************************************************
	% Ensamblado de la matriz
	%*************************************************
	A(kk:kk+1,kk:kk+1)=A(kk:kk+1,kk:kk+1)+Ak;
	%*************************************************
	% Calculo 2º miembro elemental
	%*************************************************
	ck=f*h/2*[1;1];
	%*************************************************
	% Ensamblado del segundo miembro
	%*************************************************
	c(kk:kk+1)=c(kk:kk+1)+ck;
end
%*************************************************
% Tratamiento de las condiciones de contorno
%*************************************************
%  Dirichlet en x=a:
% Bloqueo
		A(1,:)=[1,zeros(1,nodos-1)];c(1)=alfa;
        %  Dirichlet en x=b:
% Bloqueo
        A(nodos,:)=[zeros(1,nodos-1),1];c(nodos)=beta;
%*************************************************
% Resolución del sistema
%*************************************************
uh=A\c ;x=x';






