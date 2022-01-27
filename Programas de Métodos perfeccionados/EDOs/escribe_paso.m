function escribe_paso(fid,k,x,y,err)
%
% Funcion escribe_paso(fid,k,x,y,err)
%
% Funcion que imprime la evolucion de la variable independiente,
% de las variables dependientes y del error cometido (opcional)
% en los nodos de la discretizacion
% para el problema de valor inicial
%   | y'(x) = f(x,y(x))
%   | y(x=x0) = y_0
%
% Entrada:
% fid: identificador de fichero. Para salida por pantalla usar 1.
%   k: entero que indica el paso en el que estamos.
%   x: variable independiente en la EDO.
%   y: vector que contiene las variables dependientes.
% err: vector que contiene el error entre la solucion exacta y la
%      aproximada.
%
% Salida:
%   Esta rutina imprimimos por pantalla o por fichero el paso, la
%   variable independiente, las variables dependientes y, en caso
%   de que queramos, el error. Se usan formatos diferentes
%   dependiendo de la longitud de y.
%

n = length(y);

fprintf(fid,'| %5i   | %10.3e | %10.3e |',k,x,y(1));
if(n == 2)
    fprintf(fid,' %10.3e |',y(2));
end
if(nargin==5)
    fprintf(fid,' %10.3e |',err(1));
    if(n == 2)
        fprintf(fid,' %10.3e |',err(2));
    end
end

if(n>2)
    for i=2:n;
        fprintf(fid,'\n');
        fprintf(fid,'|         |            ');
        fprintf(fid,'| %10.3e |',y(i));
        if(nargin==5)
            fprintf(fid,' %10.3e |',err(i));
        end
    end
end

fprintf(fid,'\n');
fprintf(fid,'-------------------------------------');
if(n == 2 )
    fprintf(fid,'-------------');
end
if(nargin==5)
    fprintf(fid,'-------------');
    if(n == 2 )
        fprintf(fid,'-------------');
    end
end
fprintf(fid,'\n');

return
