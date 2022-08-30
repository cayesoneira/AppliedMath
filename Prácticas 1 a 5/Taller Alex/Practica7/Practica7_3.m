% PRÁCTICA 7 ejercicio 3
% Ejecutar primero la función ejer3.m
% Exporta a la ventana de comandos la MALLA y la SOLUCIÓN

whos % indica las dimensiones de las distintas matrices disponibles

% Cálculo del grado de aproximación de la solución en norma 2
error_u=u-(1-p(1,:)'.^2-p(2,:)'.^2)/4;
error2_u=norm(error_u,2)

% Cálculo del grado de aproximación de la aproximación norma infinito
errorInf_u=norm(error_u, Inf)

% Cálculo del gradiente de la solución asociada a la malla considerada 
[ux,uy]=pdegrad(p,t,u);
% obsérvese que se asigna un gradiente constante por elemento, que aproxima
% el gradiente en el baricentro del elemento

% Cálculo del gradiente correspondiente a la solución exacta por elemento 
% Se necesita primero calcular el baricentro de cada elemento
nel=size(t,2);
bar_e=[];
for k=1:nel
    bar_e=[bar_e, 1/3*[sum(p(1,t(1:3,k))), sum(p(2,t(1:3,k)))]'];
end
error2_gradu=max(norm(ux+bar_e(1,:)/2,2), norm(uy+bar_e(2,:)/2,2))
errorInf_gradu=max(norm(ux+bar_e(1,:)/2,'inf'), norm(uy+bar_e(2,:)/2,'inf'))

% Parámetro de refinamiento de malla
h=paraRef(t,p) % triángulo de área máxima

% Norma L^2
% % u_exac=@(x,y) (1-x.^2-y.^2)/4;
% % u_exacn=zeros(1,length(p(1,:)))
% % for i=1:length(p(1,:))
% %      u_exacn(i)=u_exac(p(1,i),p(2,i));
% % end
% % u_exacn=u_exacn';
% % norma_L2=normal2(u-u_exacn,t,p)
norma_L2=normal2(error_u,t,p)