function u = ef(F, m, bc, qd, varargin)
%EF  Metodo de los elementos finitos (EF) para una EDO lineal. 
%    u = ef(F, m, bc, qd, varargin) aproxima la solucion de:
%       (F.p u', v') + (F.q u', v) + (F.r u, v) = (F.f, v)
%    con condiciones de contorno
%    u(a) = bc.ua  o  -pu'(a) + bc.ha u(a) = bc.ga
%    u(b) = bc.ub  o   pu'(b) + bc.hb u(b) = bc.gb
%    sobre la malla m, con EF de grado m.deg, quadratura qd y con bloqueo
% indicado en varargin (en su defecto, condensacion).
% OBSERVACION: el termino (q u',v) queda como ejercicio 

%base de EF
switch m.deg
    case 1
        P  = @(x) [1-x, x];
        DP = @(x) [ -1, 1];
    case 2
        P  = @(x) [2*x^2-3*x+1, 4*x*(1-x), 2*x*(x-1/2)];
        DP = @(x) [      4*x-3,     4-8*x,       4*x-1];
end
%alojamiento del sistema
nnz = 3*m.nnod-2 + (m.deg==2)*2*m.nel;
A = spalloc(m.nnod, m.nnod, nnz);
b = zeros(m.nnod, 1);
%matrices elementales
for k = 1:m.nel
    Ak   = zeros(m.deg+1);
    bk   = zeros(m.deg+1, 1);
    hk   = m.xv(k+1) - m.xv(k);
    Fk   = @(x) m.xv(k) + hk*x;
    nsdk = m.nsd(k);
    for l = 1:qd.n
        xl = Fk(qd.x(l));
        Ak = Ak + qd.w(l)*( F.p(xl,nsdk)*DP(qd.x(l))'*DP(qd.x(l))/hk + ...
                            F.q(xl,nsdk)* P(qd.x(l))'*DP(qd.x(l)) + ...   
                            F.r(xl,nsdk)* P(qd.x(l))'* P(qd.x(l))*hk ); %falta (q u',v)
        bk = bk + qd.w(l)*( F.f(xl,nsdk)* P(qd.x(l))' )*hk;
    end
    %ensamblado y sin usar las matrices de ensamblado!
    n1 = m.deg*(k-1) + 1; % El primer nodo del lineal del elemento k es el k; el primer nodo del cuadrático es el 2*k-1
    nk = n1:(n1 + m.deg);
    A(nk, nk) = A(nk, nk) + Ak; %#ok<SPRIX>
    b(nk)     = b(nk)     + bk;
end
%condiciones de contorno
if isfield(bc,'ha') && isfield(bc,'ga')
    A(1,1) = A(1,1) + bc.ha;
    b(1)   = b(1)   + bc.ga;
end
if isfield(bc,'hb') && isfield(bc,'gb')
    A(end,end) = A(end,end) + bc.hb;
    b(end)     = b(end)     + bc.gb;
end
if length(varargin) == 1
    %bloqueo indicado en varargin
    [A, b] = varargin{1}(A, b, bc);
else
    
    %bloqueo por defecto: condensacion
    if isfield(bc,'ua')
        b = b(2:end) - A(2:end, 1) * bc.ua;
        A = A(2:end, 2:end);
    end
    if isfield(bc,'ub')
        b = b(1:end-1) - A(1:end-1, end) * bc.ub;
        A = A(1:end-1, 1:end-1);
    end
end
%resolucion
u = (A\b)';
if length(varargin) ~= 1
    if isfield(bc,'ua'); u = [bc.ua u]; end 
    if isfield(bc,'ub'); u = [u bc.ub]; end
end
cond = condest(A);
% disp(['El numero de condicionamiento de la matriz es ',num2str(cond)])
%2021 Francisco Pena