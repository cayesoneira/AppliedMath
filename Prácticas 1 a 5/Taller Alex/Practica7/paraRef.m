function z=paraRef(t,p)
% Función que calcula el parámetro de refinamiento de la malla como el máximo
% de las áreas de los elementos
% El cálculo se realiza sobre la malla t en el formato de pdetool
% El argumento p representa las coordenadas de los nodos

nel=length(t(1,:));
areak=ones(1,nel);
for k=1:nel
    Bk=[p(:,t(2,k))-p(:,t(1,k)), p(:,t(3,k))-p(:,t(1,k))];
    areak(k)=abs(det(Bk))*0.5;
end

z=norm(areak, Inf)