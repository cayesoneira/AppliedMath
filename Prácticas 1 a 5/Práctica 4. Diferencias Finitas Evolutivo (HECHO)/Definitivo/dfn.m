function [x, y, k] = dfn(F, bc, n, nl)
%DFIF   Metodo de diferencias finitas (iteracion funcional) para una EDO no lineal. 
%    [x,y,k] = dfif(F,bc,n,nl) aproxima la solucion con n intervalos de:
%        y'' = F.u + F.v y^F.m + F.w y'
%        condiciones de contorno mixtas (ver df.m)
%    en k iteraciones, con un iterante inicial nl.yk, un maximo de nl.maxit
%    iteraciones y tolerancia nl.tol
Fl.u = @uN; Fl.v = @vN; Fl.w = F.w;
d1 = @(y,h) (y(3:end)-y(1:end-2))/(2*h);            %aprox. primera derivada
d2 = @(y,h) (y(3:end)-2*y(2:end-1)+y(1:end-2))/h^2; %aprox. segunda derivada
for k = 1:nl.maxit
  [x, y] = df(Fl, bc, n);
  h      = x(2) - x(1);
  pc     = 2:length(x)-1;
  ydif   = norm((y-nl.yk)./max(abs(y),eps), Inf);
  Fres   = norm(d2(y,h) - F.u(x(pc)) - F.v(x(pc)).*y(pc).^F.m ...
           - F.w(x(pc)).*d1(y,h), Inf);
  if isfield(F,'lt')==0
      if Fres < nl.tol
          return
      end
  end

  if ydif < nl.tol
      return
  end
  nl.yk = y;
end
% Hice un cambio sobre el dfif original: que me avise de cuántas
% iteraciones hizo antes de fallar.
error(['El algoritmo no converge en ' num2str(nl.maxit) ' iteraciones.'])

function y = uN(x)
    y = F.u(x) + (1-F.m).* F.v(x).* nl.yk.^(F.m);
    if isfield(F,'ltY')==1
        y = F.u(x) + (1-F.m).* F.v(x).* nl.yk.^(F.m)-F.ltY(x);
    end
end


function y = vN(x)
    y = F.m.* F.v(x) .* nl.yk.^(F.m-1);
    if isfield(F,'lt')==1
        y = F.m.* F.v(x) .* nl.yk.^(F.m-1) + F.lt(x);
    end
end

end
%2021 Francisco Pena modificado por Cayetano Soneira en 2022