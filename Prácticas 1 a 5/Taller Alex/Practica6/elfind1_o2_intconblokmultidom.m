function [x,uh]=elfind1_o2_intconblokmultidom(fp,fq,ff,datos_dom,datos_a,datos_b,N,k,iopcoef,iopblo)
% [x,uh]=elfind1_o2_intconblok(fp,fq,ff,a,b,datos_a;datos_b,N,k,iopcoef,iopblo)
% elfind1_o2 calcula la soluci�n discreta del problema 
% de contorno lineal de segundo orden:
%****************************************************
% 				-(p u'(x))'+q u(x)=f  en (a,b)
% 				u(a)=alfa / p(a) u'(a)=ga / p(a) u'(a)+beta u(a)=ga 
%               u(b)=alfa /p(b) u'(b)=gb / p u'(b)+beta u(b)=gb, 
% Se suponen los coeficientes p,q,f pueden ser constantes o funciones
%   las variables alfa y beta son constantes
% Se considera que el coeficiente p puede ser distinto en cada subdominio.
%****************************************************
% utilizando elementos finitos de Lagrange de grado k=1.
% N+1 n�mero de elementos de la malla.
% datos_dom es un vector de dimensi�n n� de subdominios+2
%   datos_dom(1)=n�mero de subdominios
%   datos_dom(i)=extremo izquierdo del subdominio i-1, 2<=i<=n� de
%   subdominios+1
%   datos_dom(n� de subdominios + 2)=extremo derecho del intervalo total.
%k: orden de los polinomios de interpolaci�n;
%   k puede ser 1 � 2
% N : n�mero de puntos de discretizaci�n interiores al dominio.
% iopcoef: par�metro de elecci�n del m�todo de integraci�n si k=1
%           =1 si los coeficientes son constantes
%           =2 si se usa integraci�n de Poncelet
%           =3 si se usa la f�rmula del trapecio
%datos_a: permiten elegir la condici�n de borde en a
%           =[1,alfa] si es Dirichlet;
%           =[2,ga] si es Neumann
%           =[3,alfa,ga] si es Robin
%           Una definici�n an�loga para datos_b
%iopblo: opci�n de bloqueo
%       =1 se reemplaza la correspondiente ecuaci�n del sistema por
%               u(a)=alfa /u(b)=alfa
%       =2 se reemplaza la correspondiente ecuaci�n del sistema por
%               1.e30*u_1+a(1,j)*u(j)=alfa /a(n,j)*u(j)+1.e30*u_nodos=alfa
%       =3 se utiliza una t�cnica de condensaci�n
%****************************************************
% Construcci�n de la malla 1-dimensional.
%****************************************************
global h nodos
% Inicializaci�n de los valores de h, elementos, y longitud para cada
% subdominio
ndom=datos_dom(1);
h=zeros(1,ndom); % tama�o de la discretizacion;
nel1=zeros(1,ndom);
long=zeros(1,ndom);
nver=0;
x=[];

if ndom ==1
    % c�lculo de los extremos del intervalo
    a=datos_dom(2);b=datos_dom(3);
    h(1)=(b-a)/(N+1); % tama�o de la discretizacion
    nel1(1)=N+1; % numero de elementos: T_0,T_1,...,T_N
    nver=N+2; % numero de v�rtices: a_0,a_1,...,a_N+1
    i=1:nver; x=a+(i-1)*h(1);%coordenadas de los v�rtices
else
    for j=1: ndom
    a=datos_dom(2);b=datos_dom(ndom+2);
    h(j)=(b-a)/(N+1); % tama�o inicial que corresponder�a para una discretizacion regular
    long(j)=datos_dom(j+2)-datos_dom(j+1);
    nel1(j)=round(long(j)*(N+1)/(b-a));
    h(j)=long(j)/nel1(j); % tama�o de elemento en el dominio j para que la malla sea conforme.
    i=1:nel1(j); x=[x, datos_dom(j+1)+(i-1)*h(j)];
    nver=nver+nel1(j);
    end
    x=[x, b];
    nver=nver+1;
   
end
    nel=sum(nel1); 
    
   
M=k*nel+1-2;
nodos=M+2;%n�mero total de nodos o grados de libertad

%***********************************************
% Inicializaci�n de la matriz A
%***********************************************
if k==1
    nonulos=3*nodos-2;
else
    nonulos=3*nodos-2+2*nel;
    x_nodos=[x(1)]; % inicializaci�n vector de coordenadas de los nodos
end

A=spalloc(nodos,nodos,nonulos);
% Inicializaci�n de c
%************************************************
c=zeros(nodos,1);
%*************************************************
% Bucle en elementos
%*************************************************
nelj=0;
for j=1:ndom
        
for kk=nelj+1:nelj+nel1(j)
    %*************************************************
    % C�lculo de los coeficientes en el elemento kk
    % C�lculo submatrices elementales
	% C�lculo 2� miembro elemental
    %*************************************************
    if k==1
        if iopcoef==1
    Kk=fp/h(j)*[1, -1;-1, 1];Mk=(fq*h(j))/6*[2 1;1 2];
    ck=ff*h(j)/2*[1;1];
    elseif iopcoef==2
    punto_medio=(x(kk)+x(kk+1))/2;
    p=feval(fp,punto_medio,j);q=feval(fq,punto_medio);
    Kk=p/h(j)*[1, -1;-1, 1];Mk=(q*h(j))/4*[1 1;1 1];
    f=feval(ff,punto_medio);
    ck=f*h(j)/2*[1;1];
    elseif iopcoef==3
    p1=feval(fp,x(kk),j);p2=feval(fp,x(kk+1),j);
    q1=feval(fq,x(kk));q2=feval(fq,x(kk+1));
    Kk=(p1+p2)/2/h(j)*[1, -1;-1, 1];Mk=h(j)/2*[q1 0;0 q2]; 
    f1=feval(ff,x(kk));f2=feval(ff,x(kk+1));
    ck=h(j)/2*[f1;f2];
        end
    elseif k==2
        punto_medio=(x(kk)+x(kk+1))/2;
        x_nodos=[x_nodos, punto_medio,x(kk+1)];
        Mk=1/(6*h(j))*([9 -12 3;-12 16 -4; 3 -4 1]*feval(fp,x(kk),j)+...
            4*[1 0 -1; 0 0 0; -1 0 1]*feval(fp,punto_medio,j)+...
            [1 -4 3; -4 16 -12; 3 -12 9]*feval(fp,x(kk+1),j));
        Kk=h(j)/6*([feval(fq,x(kk)) 0 0;0 4*feval(fq,punto_medio) 0; 0 0 feval(fq,x(kk+1))]);
        if iopcoef==1
         ck=ff*h(j)/6*[1;4;1];  
        else
            ck=h(j)/6*([feval(ff,x(kk)) ; 4*feval(ff,punto_medio); feval(ff,x(kk+1))]);
        end
    end
    %*************************************************
	% Calculo matriz elemental
	%*************************************************
   Ak=Mk+Kk;
	%*************************************************
	% Ensamblado de la matriz
	%*************************************************
    k1=k*(kk-1)+1;
	A(k1:k1+k,k1:k1+k)=A(k1:k1+k,k1:k1+k)+Ak;
	%*************************************************
	% Ensamblado del segundo miembro
	%*************************************************
	c(k1:k1+k)=c(k1:k1+k)+ck;
end
nelj=nelj+nel1(j);
end

%*************************************************
% Tratamiento de las condiciones de contorno
%*************************************************
if datos_a(1)==1 %Dirichlet en x=a:
        if(iopblo==1)
            A(1,:)=[1,zeros(1,nodos-1)];c(1)=datos_a(2);
        elseif(iopblo==2)
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
% Resoluci�n del sistema
%*************************************************
uh=A\c ;
if k==1
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







