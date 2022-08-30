function [x,uh]=elfin1dns(fp,fdp,fr,fq,ff,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo)
% [x,uh]=elfin1d(fp,fq,ff,a,b,datos_a;datos_b,N,k,iopcoef,iopblo)
% elfind1_o2 calcula la solución discreta del problema 
% de contorno lineal de segundo orden:
%****************************************************
% 				-p u''(x)+ru'(x)+q u(x)=f  en (a,b)
% 				u(a)=alfa / p(a) u'(a)=ga / p(a) u'(a)+alfa u(a)=ga 
%               u(b)=alfa /p(b) u'(b)=gb / p(b) u'(b)+beta u(b)=gb, 
% Se suponen los coeficientes p,q,f pueden ser constantes o funciones
% Los argumentos alfa y beta son constantes
%****************************************************
% utilizando elementos finitos de Lagrange de grado k=1.
% nel número de elementos de la malla.
%k: orden de los polinomios de interpolación;
%   l puede ser 1 ó 2
% iopcoef: parámetro de elección del método de integración si l=1
%           =1 si los coeficientes son constantes
%           =2 si se usa integración de Poncelet
%           =3 si se usa la fórmula del trapecio
% Si l=2 se utiliza la fórmula de Simpson
%datos_a: permiten elegir la condición de borde en a
%           =[1,alfa] si es Dirichlet;
%           =[2,ga] si es Neumann
%           =[3,alfa,ga] si es Robin
%           Una definición análoga para datos_b
%iopblo: opción de bloqueo
%       =1 se reemplaza la correspondiente ecuación del sistema por
%               u(a)=alfa /u(b)=alfa
%       =2 se reemplaza la correspondiente ecuación del sistema por
%               1.e30*u_1+a(1,j)*u(j)=alfa /a(n,j)*u(j)+1.e30*u_nodos=alfa
%       =3 se utiliza una técnica de condensación
%****************************************************
% Construcción de la malla 1-dimensional.
%****************************************************
global h nodos
nver=nel+1; % numero de vértices
h=(b-a)/nel; % tamaño de la discretizacion
nodos=l*nel+1;%número total de nodos o grados de libertad
i=1:nver; x=a+(i-1)*h;%coordenadas de los vértices
%***********************************************
% Inicialización de la matriz A
%***********************************************
if l==1
    nonulos=3*nodos-2;
else
    nonulos=3*nodos-2+2*nel;
    x_nodos=[x(1)]; % inicialización vector de coordenadas de los nodos
end

A=spalloc(nodos,nodos,nonulos);
% Inicialización de c
%************************************************
c=zeros(nodos,1);
%*************************************************
% Bucle en elementos
%*************************************************
for k=1:nel
    %*************************************************
    % Cálculo de los coeficientes en el elemento
    % Calculo submatrices elementales
	% Calculo 2º miembro elemental
    %*************************************************
%     if l==1
        if iopcoef==1
    Kk=fp/h*[1, -1;-1, 1];Mk=(fq*h)/6*[2 1;1 2];Rk=fr/2*[-1 1;-1 1];
    ck=ff*h/2*[1;1];
    elseif iopcoef==2
    punto_medio=(x(k)+x(k+1))/2;
    p=feval(fp,punto_medio);q=feval(fq,punto_medio);dpr=feval(fdp,punto_medio)+feval(fr,punto_medio);
    Kk=p/h*[1, -1;-1, 1];Mk=(q*h)/4*[1 1;1 1];Rk=dpr/2*[-1 1;-1 1];
    f=feval(ff,punto_medio);
    ck=f*h/2*[1;1];
    elseif iopcoef==3
    p1=feval(fp,x(k));p2=feval(fp,x(k+1));
    q1=feval(fq,x(k));q2=feval(fq,x(k+1));
    dpr1=feval(fdp,x(k))+feval(fr,x(k));
    dpr2=feval(fdp,x(k+1))+feval(fr,x(k+1));
    Kk=(p1+p2)/2/h*[1, -1;-1, 1];Mk=h/2*[q1 0;0 q2]; Rk=1/2*[-dpr1 dpr1;-dpr2 dpr2];
    f1=feval(ff,x(k));f2=feval(ff,x(k+1));
    ck=h/2*[f1;f2];
        end
   % elseif l==2
   %     punto_medio=(x(k)+x(k+1))/2;
   %    x_nodos=[x_nodos, punto_medio,x(k+1)];
%         Mk=1/(6*h)*([9 -12 3;-12 16 -4; 3 -4 1]*feval(fp,x(k))+4*[1 0 -1; 0 0 0; -1 0 1]*feval(fp,punto_medio)+[1 -4 3; -4 16 -12; 3 -12 9]*feval(fp,x(k+1)));
%         Kk=h/6*([feval(fq,x(k)) 0 0;0 4*feval(fq,punto_medio) 0; 0 0 feval(fq,x(k+1))]);
%         if iopcoef==1
%          ck=ff*h/6*[1;4;1];  
%         else
%             ck=h/6*([feval(ff,x(k)) ; 4*feval(ff,punto_medio); feval(ff,x(k+1))]);
%         end
%     end
    %*************************************************
	% Calculo matriz elemental
	%*************************************************
   Ak=Mk+Kk+Rk;
	%*************************************************
	% Ensamblado de la matriz
	%*************************************************
    k1=l*(k-1)+1;
	A(k1:k1+l,k1:k1+l)=A(k1:k1+l,k1:k1+l)+Ak;
	%*************************************************
	% Ensamblado del segundo miembro
	%*************************************************
	c(k1:k1+l)=c(k1:k1+l)+ck;
end
%*************************************************
% Tratamiento de las condiciones de contorno
%*************************************************
if datos_a(1)==1 %Dirichlet en x=a:
        if(iopblo==1)
            A(1,:)=[1,zeros(1,nodos-1)];c(1)=datos_a(2);
        elseif(iopblo==2)
            % c(1)=1.e30*A(1,1)*datos_a(2);
            % A(1,1)=1.e30*A(1,1);
		A(1,1)=1.e30;c(1)=1.e30*datos_a(2);
        elseif(iopblo==3)
           c= c(2:nodos)-A(2:nodos,1)*datos_a(2);
           A=A(2:nodos,2:nodos);
           nodos=nodos-1;
        end
elseif datos_a(1)==2 %Neumann en x=a       
        c(1)=c(1)-datos_a(2);
elseif datos_a(1)==3 %Robin en x=a
        A(1,1)=A(1,1)-datos_a(2);
        c(1)=c(1)-datos_a(3);
end
if datos_b(1)==1 %Dirichlet en x=b:
    if(iopblo==1)
            A(nodos,:)=[zeros(1,nodos-1),1];c(nodos)=datos_b(2);
        elseif(iopblo==2)
		A(nodos,nodos)=1.e30;c(nodos)=1.e30*datos_b(2);
        elseif(iopblo==3)
           c= c(1:nodos-1)-A(1:nodos-1,nodos)*datos_b(2);
           A=A(1:nodos-1,1:nodos-1);
           nodos=nodos-1;
    end
elseif datos_b(1)==2 %Neumann en x=b       
        c(nodos)=c(nodos)+datos_b(2);
elseif datos_b(1)==3 %Robin en x=b
        A(nodos,nodos)=A(nodos,nodos)+datos_b(2);
        c(nodos)=c(nodos)+datos_b(3);
end
   
%*************************************************
% Resolución del sistema
%*************************************************
uh=A\c ;
if l==1
    x=x';
else
    x=x_nodos';
end

if(iopblo==3)&& datos_a(1)==1
    uh=[datos_a(2);uh];
    nodos=nodos+1;
end
if(iopblo==3)&& datos_b(1)==1
    uh=[uh;datos_b(2)];
    nodos=nodos+1;
end







