function escribe_cabecera(fid,x,y,err)
%
% Funcion escribe_cabecera(fid,x,y,err)
%
% Funcion que imprime un encabezado de cara a imprimir 
% la evolucion de la variable independiente, 
% de las variables dependientes y del error cometido (opcional) 
% en los nodos de la discretizacion
% para el problema de valor inicial
%   | y'(x) = f(x,y(x))
%   | y(x=x0) = y_0
% Se trata de una rutina que se debe llamar antes de escribe_paso.m
%
% Entrada:
% fid: identificador de fichero. Para salida por pantalla usar 1.
%   x: variable independiente en la EDO (no se usa en esta funcion).
%   y: vector que contiene las variables dependientes (su utilidad se 
%      limita a conocer el numero de componentes de las que dispone el 
%      el vector y).
% err: vector que contiene el error entre la solucion exacta y la 
%      aproximada (no se usa en esta funcion).
%
% Salida: 
%   Esta rutina tan solo proporciona un encabezado cuyo formato dependera
%   del numero de componentes en y. La salida es o bien por fichero o por 
%   pantalla.
%   

n = length(y);

fprintf(fid,'-------------------------------------');
if(n == 2 )
    fprintf(fid,'-------------');
end
if(nargin==4)
    fprintf(fid,'-------------');
    if(n == 2 )
        fprintf(fid,'-------------');
    end
end
fprintf(fid,'\n');
fprintf(fid,'|  Nodo   |     x      |');
if(n > 2 || n == 1)
    fprintf(fid,       '      y     |');
elseif (n == 2)
    fprintf(fid,       '     y_1    |     y_2    |');
end
if(nargin==4)
    if(n > 2 || n == 1)
        fprintf(fid,   '    err     |');
    elseif (n == 2)
        fprintf(fid,   '    err_1   |    err_2   |');
    end
end
fprintf(fid,'\n');
fprintf(fid,'-------------------------------------');
if(n == 2 )
    fprintf(fid,'-------------');
end
if(nargin==4)
    fprintf(fid,'-------------');
    if(n == 2 )
        fprintf(fid,'-------------');
    end
end

fprintf(fid,'\n');

return 
