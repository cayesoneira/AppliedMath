function [x, y] = df(F, bc, n)
%DF   Metodo de diferencias finitas para una EDO lineal con condiciones mixtas. 
%    [x,y] = df(F,bc,n) aproxima la solucion en n pasos (en las abscisas x) de:
%        y'' = F.u + F.v y + F.w y' 
%        y(bc.a) = bc.ya si bc.blk_a, si no -y'(bc.a) + bc.ha y(bc.a) = bc.ya 
%        y(bc.b) = bc.yb si bc.blk_b, si no  y'(bc.b) + bc.hb y(bc.b) = bc.yb
h  = (bc.b - bc.a) / n;
x  = bc.a + h*(0:n);
ux = F.u(x); vx = F.v(x); wx = F.w(x);
%matriz:      filas cols         valor           size 
A =     sparse(2:n, 2:n,    2 + h^2*vx(2:n),   n+1, n+1); %diag
A = A + sparse(2:n, 1:n-1, -1 - h  *wx(2:n)/2, n+1, n+1); %subdiag
A = A + sparse(2:n, 3:n+1, -1 + h  *wx(2:n)/2, n+1, n+1); %superdiag
%segundo miembro
b(2:n) = -h^2*ux(2:n);
%condiciones de contorno
if bc.blk_a
    A(1,1) = 1;
    b(1)   = bc.ya;
else
    A(1,1:2) = [1 + h*bc.ha + h^2*( vx(1) + bc.ha*wx(1))/2, -1];
    b(1)     = h*bc.ya + h^2*(-ux(1) + bc.ya*wx(1))/2;
end
if bc.blk_b
    A(n+1,n+1) = 1;
    b(n+1)     = bc.yb;
else
    A(n+1,n:n+1) = [-1, 1 + h*bc.hb + h^2*( vx(n+1) - bc.hb*wx(n+1))/2];
    b(n+1)       = h*bc.yb - h^2*(ux(n+1) + bc.yb*wx(n+1))/2;
end
%solucion
y = (A\b')';
%2021 Francisco Pena