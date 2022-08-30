function [x, t, Y] = dfnle_esquema(F, bi, nx, nt, nl, algor)
% DFIF   Metodo de diferencias finitas para una EDO no lineal evolutiva. 
%    [x,t,y] = dfnle(f,bi,nx,nt,nl) aproxima la solucion con nx,nt intervalos en (x,t) de:
%        F.l dy/dt - d2y/dx2 + F.u + F.v y^F.m + F.w dy/dx = 0 en [a,b] X [t0,tf]
%        y(a,t) = bi.ya(t) si bi.blk_a, si no -dy/dx(a,t) = bi.ya(T) - bi.ha Y(a,t)
%        y(b,t) = bi.yb(t) si bi.blk_b, si no  dy/dx(b,t) = bi.yb(T) - bi.hb Y(b,t)
%        y(x,bi.t0) = bi.y0 
%    en K iteraciones del algoritmo algor, con un iterante inicial nl.yk, 
%    un maximo de nl.maxit iteraciones y tolerancia nl.tol
bc =  bi; %copia a, b, blk_a, blk_b, ha, hb
h  = (bi.b  - bi.a) /nx;
dt = (bi.tf - bi.t0)/nt;
x  =  bi.a  +  h*(0:nx);
t  =  bi.t0 + dt*(0:nt);
Y  = zeros(nt+1, nx+1);
Y(1,:) = bi.y0(x);
nl.yk = Y(1,:);
for j = 2:length(t)
    bc.ya = bi.ya(t(j));
    bc.yb = bi.yb(t(j));
    if F.m == 1
        Fs = struct('u', @(s)  F.u(s,t(j)) - F.l(s,t(j)) .* Y(j-1,:) / dt, ...
                    'v', @(s)  F.v(s,t(j)) + F.l(s,t(j)) ./ dt, ...
                    'w', @(s)  F.w(s,t(j)));

        [~, Y(j,:)] = algor(Fs, bc, nx);
        
%         tiledlayout('flow')
%         nexttile
%         plot(x,Y(j,:))

    else
        % JUSTO AQUÍ ADELANTE FALLA PORQUE Y(j-1,:) NO SIRVE. LO CAMBIO POR
        % ONES(sze(x))
        Fs = struct('ltY', @(s)  F.l(s,t(j)).* Y(j-1,:) ./bc.dt, ...
                    'lt', @(s) F.l(s,t(j)) ./ bc.dt, ...
                    'u', @(s)  F.u(s,t(j)), ...
                    'v', @(s)  F.v(s,t(j)), ...
                    'w', @(s)  F.w(s,t(j)), ...
                    'm', F.m);

        nl.yk = Y(j-1,:);
        [~, Y(j,:)] = algor(Fs, bc, nx, nl);
    end
end
%2021 Francisco Pena