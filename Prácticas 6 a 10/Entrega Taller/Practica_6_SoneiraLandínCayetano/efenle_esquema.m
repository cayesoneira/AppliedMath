function [t, U] = efenle(F, m, bi, qd, nl, algor, varargin)
%EFIF  Metodo de los elementos finitos para una EDO no lineal evolutiva.
%    [T, U] = efnle(F,M,BI,NT,QD,NL,ALGOR,VARARGIN) aproxima la solucion con NT intervalos en T de:
%        (F.L dU/dt, V) + (F.P dU/dx, dV/dx) + (F.Q dU/dx, V) + (F.R U^F.M, V) = (F.F, V)
%        U(A,T) = BI.UA(T) si BI.BLK_A, si no -dU/dX(A,T) = BI.UA(T) - BC.HA U(A,T)
%        U(B,T) = BI.UB(T) si BI.BLK_B, si no  dU/dX(B,T) = BI.UB(T) - BC.HB U(B,T)
%        U(X,BI.T0) = BI.U0
%    con condiciones de contorno BC, malla M, cuadratura QD y bloqueo
%    opcional VARARGIN (ver EF.M). Iteracion funcional aproxima la solucion en en K iteraciones
%    del algoritmo ALGOR, con un iterante inicial NL.UK, un maximo de NL.MAXIT iteraciones
%    y tolerancia NL.TOL

bc = bi; %copia a, b, ha, hb
nt=bi.nt;
dt = (bi.tf - bi.t0)/nt;
t  = bi.t0 + dt*(0:nt);
U = zeros(nt+1, m.nnod);
U(1,:) = bi.u0;
for j = 2:length(t)
    if isfield(bc,'ua'), bc.ua = bi.ua(t(j)); end
    if isfield(bc,'ub'), bc.ub = bi.ub(t(j)); end
    if isfield(bc,'ga'), bc.ga = bi.ga(t(j)); end
    if isfield(bc,'gb'), bc.gb = bi.gb(t(j)); end

    if F.m == 1
        Fs = struct('p', @(x,nsd) F.p(x,t(j)),...
            'q', @(x,nsd) F.q(x,t(j)),...
            'r', @(x,nsd) F.r(x,t(j))+ interp1(m.xn,U(j-1,:),x)) * F.lt(x,t(j))/d,...
            'f', @(x,nsd) F.f(x,t(j)) + F.l(x,t(j))*interp1(m.xn,U(j-1,:),x)/dt;
        U(j,:) = ef(Fs, m, bc, qd, varargin{:});
    else
        Fs = struct('lt', @(x,nsd) F.l(x,t(j)),...
            'p', @(x,nsd) F.p(x,t(j)),...
            'q', @(x,nsd) F.q(x,t(j)),...
            'r', @(x,nsd) F.r(x,t(j)),...
            'f', @(x,nsd) F.f(x,t(j)),...
            'm',  F.m);
        nl.uk = U(j-1,:);
        U(j,:) = algor(Fs, m, bc, qd, nl, varargin{:});
    end
end
%2021 Francisco Pena