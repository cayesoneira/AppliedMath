function [t,x,h]=diffinc(u,v,w,a0,b0,alfa,beta,npas,iopc,gamma,delta)
% [t,x,h]=difinc(u,v,w,a0,b0,alfa,beta,num_pas,iopc,gamma,delta)
%       diffinc calcula la soluci�n discreta del problema de contorno
%       lineal de segundo orden:
%   x''(t)=u(t)+v(t)x+w(t)x'
% donde las condiciones de contorno dependen del valor de iopc
%   iopc=1 ........x(a)=alfa, x(b)=beta
%   iopc=2 ........x'(a)+gamma x(a)=alfa, x(b)=beta
%   iopc=3 ........x(a)=alfa,x'(b)+delta x(b)=beta 
%   iopc=4 ........x'(a)+gamma x(a)=alfa, x'(b)+delta x(b)=beta
% ARGUMENTOS DE ENTRADA:
% u,v,w : strings que indican los coeficientes de la ecuaci�n.
% a0,b0 : extremos del intervalo sobre el que se integra la ecuaci�n.
% alfa,beta : Valor de la condici�n de contorno en ambos extremos.
% iopc : argumento opcional, si no existe se resolver� el problema
%        Dirichlet en ambos extremos.
%        Si iopc existe jugar� el papael indicado m�s arriba.
% gamma, delta: argumentos opcionales que deben existir si iopc existe
%              y es distinto de 1. Indicar�n los coeficientes 
%              de las condiciones Robin.
% npas : n�mero de pasos, correspondiente a la variable de las notas
%           num_pas
% ARGUMENTOS DE SALIDA:
% t : vector fila que contiene los puntos de discretizaci�n de la variable
%       independiente.
% x : vector fila que contiene el valor aproximado en los puntos de discretizaci�n 
%       para  la variable dependiente.
% h : argumanto de salida opcional, que indica el par�metro de
%       discretizaci�n h.
%__________________________________________________________________________
% An�lisis argumentos de entrada
switch nargin
    case 8
        iopc=1;
    case 9     
        if (iopc~=1)
            error('Argumentos de entrada incorrectos')
        end
    case 11 
        if any(eq(iopc,[1:4]))==0
            error('Argumento de entrada iopc incorrecto')
        end
    otherwise
        error('Argumentos de entrada incorrectos')
end

%__________________________________________________________________________
% An�lisis argumentos de salida

if all(ne(nargout,[2,3]))
            error('N�mero de argumentos de salida incorrectos')
end

    
%__________________________________________________________________________
% Definici�n del par�metro de discretizaci�n.
%   h : tama�o del paso de discretizaci�n
h=(b0-a0)/npas;
%__________________________________________________________________________
% Obtenci�n del vector de variable independiente discretizado
%   i : vector de �ndices de dimensi�n npas+1
%   t : puntos de discretizaci�n
i=[0:npas];t=a0+i*h;
%__________________________________________________________________________
% C�lculo de los vectores que definen la matriz del sistema
% d : es el vector relacionado con la diagonal principal
% a : es el vector relacionado con la diagonal inferior
% c : es el vector relacionado con la diagonal superior
% b : es el vector relacionado con el segundo miembro 
% Su dimensi�n es npas-1
% a=-1*ones(1,npas+1)-h*feval(w,t)./2;
% d=2+h^2*feval(v,t);c=-a-2;
% b=-h^2*feval(u,t);
a=-1*ones(1,npas+1)-h*w(t)./2;% Opci�n para evaluar funciones de forma m�s directa.
d=2+h^2*v(t);c=-a-2;
b=-h^2*u(t);
%__________________________________________________________________________
%%C�lculo de la matriz del sistema
if iopc==1 
	matriz=sparse(1:npas-1,1:npas-1,d(2:npas));
	matriz=matriz+sparse(2:npas-1,1:npas-2,a(3:npas),npas-1,npas-1);
	matriz=matriz+sparse(1:npas-2,2:npas-1,c(2:npas-1),npas-1,npas-1);
elseif iopc==2
   matriz=sparse(2:npas,2:npas,d(2:npas),npas,npas);
   matriz=matriz+sparse(2:npas,1:npas-1,a(2:npas),npas,npas);
   matriz=matriz+sparse(2:npas-1,3:npas,c(2:npas-1),npas,npas);
   matriz(1,1)=1-gamma*h+(d(1)-2)/2-(c(1)+1)*h*gamma;
   matriz(1,2)=-1;
elseif iopc==3
   matriz=sparse(1:npas-1,1:npas-1,d(2:npas),npas,npas);
   matriz=matriz+sparse(2:npas-1,1:npas-2,a(3:npas),npas,npas);
   matriz=matriz+sparse(1:npas-1,2:npas,c(2:npas),npas,npas);
   matriz(npas,npas)=1+delta*h+(d(npas+1)-2)/2-(c(npas+1)+1)*h*delta;
   matriz(npas,npas-1)=-1;
elseif iopc==4
   matriz=sparse(2:npas,2:npas,d(2:npas),npas+1,npas+1);
   matriz=matriz+sparse(2:npas,1:npas-1,a(2:npas),npas+1,npas+1);
   matriz=matriz+sparse(2:npas,3:npas+1,c(2:npas),npas+1,npas+1);
   matriz(1,1)=1-gamma*h+(d(1)-2)/2-(c(1)+1)*h*gamma;
   matriz(1,2)=-1;
   matriz(npas+1,npas+1)=1+delta*h+(d(npas+1)-2)/2-(c(npas+1)+1)*h*delta;
   matriz(npas+1,npas)=-1;
end
 %__________________________________________________________________________
%C�lculo del segundo miembro
bmod=b(2:npas);
if(iopc==1)||(iopc==3)
   bmod(1)=bmod(1)-a(2)*alfa;
else
   bmod=[-alfa*h+b(1)/2-(c(1)+1)*h*alfa bmod];
end
if(iopc==1)
   bmod(npas-1)=bmod(npas-1)-c(npas)*beta;
elseif(iopc==2)
   bmod(npas)=bmod(npas)-c(npas)*beta;
else
   bmod=[bmod beta*h+b(npas+1)/2-(c(npas+1)+1)*h*beta];
end
   
%__________________________________________________________________________
%C�lculo de la soluci�n
x=matriz \ bmod';
if (iopc==1)
   x=[alfa,x',beta];
elseif(iopc==2)
   x=[x',beta];
elseif(iopc==3)
   x=[alfa x'];
else
   x=x';
end






