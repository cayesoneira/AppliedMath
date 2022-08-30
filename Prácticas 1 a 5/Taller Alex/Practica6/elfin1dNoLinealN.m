function [x,phi]=elfin1dNoLinealN(fp,dfp,fr,fq,ff,p,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo,tol,maxit)
global phi0 p xnodos Gamma

% Problema de contorno no lineal y no simétrico:
%   -(p(x)u'(x))'+r(x)u'(x)+q(x) u^p(x)=f 

% Para p(x) constante, la ecuación se puede reescribir como
%    u''(x)=f/p+q/p*u^p(x)+r/p*u'(x)
% El operador de Newton es 
%   q''(x)=(f/p(x)+(1-p)*q/p(x)phi0^p)+(p*q/p(x)phi0^(p-1))*q(x)+r/p(x)*q'(x)
% Para resolver en elfin1dns
%   -p(x)q''(x)+ (p*q(x)*phi0^(p-1))*q(x) + r(x)*q'(x) =-(f+(1-p)*q*phi0^p)

nver=nel+1; %nº de vértices
h=(b-a)/nel; %paso de discretización
Nodos=l*nel+1; %nº de grados de libertad o nodos
i=1:nver;x=a+(i-1)*h;
if l==1
    xnodos=x;
else
    xnodos=[];
    xmed=(x(1:nver-1)+x(2:nver))/2;
    for i=1:nver-1
        xnodos=[xnodos,x(i),xmed(i)];
    end
    xnodos=[xnodos,x(nver)];
end

phi0=ones(1,Nodos); % iterante inicial

%Coeficiente en x(t) de la ecuación linealizada
fpN=@pN;
fqN=@qN;

for i=1:maxit
    %Resolución del modelo linealizado
    [x,phi]=elfin1dns(fp,dfp,fr,fqN,fpN,a,b,datos_a,datos_b,nel,l,iopcoef,iopblo); 
    phi2=phi0+(phi0==0)*eps;
    error_rel=norm((phi-phi0)/phi2,Inf);
    
    if error_rel<tol
        fprintf(1,'Se ha alcanzado la convergencia en %5.0f iteraciones \n',i)
        return
    else
        %Actualización del coeficiente del modelo linealizado
        phi0=phi;
        global phi0
    end
end

if i==maxit
    disp('No hay convergencia')
    disp(['Nº iteraciones= ',num2str(i), ' error entre iterantes = ', num2str(error_rel)])
end
