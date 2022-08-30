function z=normal2(u,t,p,opint)
% Función que calcula la norma L^2 de una función u definida en cada nodo
% de la malla sobre el dominio definido por la malla t en el formato de
% pdetool.
% La función u debe estar definida sobre los nodos de la malla.
% La matriz p representa las coordenadas de los nodos.
% La matriz t representa la información sobre los elementos: sus vértices y
% el dominio al que pertenecen.
% opint es la opción de integración.
% opint==1, si se utiliza la aproximación numérica del valor medio
% opint==2, si se calcula la integral exacta sobre cada elemento
% Si el número de argumentos es 3, se entiende que opint=1

if nargin==3
    opint=1;
end
if opint==2
    syms x y
    PRef=[1-x-y, x, y];
    integrando=PRef'*PRef;
    Integ1=int(integrando, 'y', 0, 1-x);
    Integ2=int(Integ1,'x',0,1);
    IntegralRef=eval(Integ2);
end
nel=length(t(1,:));
z=0;
for k=1:nel
    uk=u(t(1:3,k));
    Bk=[p(:,t(2,k))-p(:,t(1,k)), p(:,t(3,k))-p(:,t(1,k))];
    if opint==1
        areak=abs(det(Bk))*0.5;
        Integralk=sum(uk.^2)/3*areak;
    else
        Integralk=abs(det(Bk))*uk'*IntegralRef*uk;
    end
    z=z+Integralk;
end
z=sqrt(z);
end