function [x,phi,i]=diffincNoLinealIF_evol(l,u,v,w,p,a,b,alfa,beta,num_pas,iopc,gamma,delta,tol,maxit)
% diffinc calcula la solución discreta del problema de contorno  lineal de
% segundo orden:
% x''=u(t)+v(t)x^p+w(t)x'
% donde las condiciones de contorno dependen del valor de iopc
%   iopc=1 ........x(a)=alfa, x(b)=beta
%   iopc=2 ........x'(a)+gamma x(a)=alfa, x(b)=beta
%   iopc=3 ........x(a)=alfa,x'(b)+delta x(b)=beta 
%   iopc=4 ........x'(a)+gamma x(a)=alfa, x'(b)+delta x(b)=beta
% ARGUMENTOS DE ENTRADA:
% u,v,w : strings que indican los coeficientes de la ecuación.
% p : exponente de la no linealidad de x(t)
% a,b : extremos del intervalo sobre el que se integra la ecuación.
% alfa,beta : Valor de la condición de contorno en ambos extremos.
% iopc : argumento opcional, si no existe se resolverá el problema
%        Dirichlet en ambos extremos.
%        Si iopc existe jugará el papael indicado más arriba.
% gamma, delta: argumentos opcionales que deben existir si iopc existe
%              y es distinto de 1. Indicarán los coeficientes 
%              de las condiciones Robin.
% npas : número de pasos
% tol : parámetro de tolerancia para la convergencia del método
% maxit : número máximo de itereaciones 

% Itereante inicial
phi0=ones(1,num_pas+1);      

% Coeficiente en x(t) de la ecuación linealizada
vIF=@(t) feval(v,t).*phi0.^(p-1)+feval(l,t);

for i=1:maxit
    % Resolución del modelo linealizado
    [x,phi]=diffinc(u,vIF,w,a,b,alfa,beta,num_pas,iopc,gamma,delta);
    phi2=phi0+(phi==0)*eps;
    error_rel=norm((phi-phi0)./phi2,Inf);
   
    if error_rel<tol
        %fprintf(1,'Se ha alcanzado la convergencia en %5.0f iteraciones \n',i)
        return
    else 
        % Actualizanción del coefeiciente del modelo linealizado
        phi0=phi;
        vIF=@(t) feval(v,t).*phi0.^(p-1)+feval(l,t);
    end 
end

