function [u, k] = efn(F, m, bc, qd, nl, varargin)
%EFIF  Metodo de los elementos finitos (Newton) para una EDO no lineal.
%    u = efn(m,bc,qd,nl,varargin) aproxima la solucion de:
%       (F.p u', v') + (F.q u', v) + (F.r u^F.m, v) = (F.f, v)
%    con condiciones de contorno bc, malla m, cuadratura qd y bloqueo
%    opcional varargin (ver ef.m). Newton usa el iterante inicial
%    nl.uk, numero maximo de iteraciones nl.maxit y tolerancia nl.tol.

Fl   = F; %copia p, q
Fl.r = @rN; 
Fl.f = @fN;
for k = 1:nl.maxit
  u = ef(Fl, m, bc, qd, varargin{:});
  % varargin{:} es la forma de tomar "todo varargin". Como es una celda...

  % Se podría añadir al criterio de parada el tamaño del residuo de F. No
  % se hizo por simplicidad.

  if norm((u-nl.uk)./max(abs(u),eps), Inf) < nl.tol
    return
  end
  nl.uk = u;
end
error('El algoritmo no converge.')

function y = fN(x, nsd)
%FN   Coeficiente f(x) para la EDO linealizada de Newton.
ukx = interp1(m.xn, nl.uk, x);
y   = F.f(x,nsd) + (F.m-1)*F.r(x,nsd)*ukx^F.m;
end

function y = rN(x,nsd)
%RN   Coeficiente r(x) para la EDO linealizada de Newton.
ukx = interp1(m.xn, nl.uk, x);
y   = F.m * F.r(x,nsd) * ukx^(F.m-1);
if isfield(F,'lt'), y = y + F.lt(x,nsd); end
end
end
%2021 Francisco Pena
