function [x, y] = dfd(F, bc, n)
%DFD   Metodo de diferencias finitas para una EDO lineal con condiciones Dirichlet. 
%    [x,y] = dfd(F,bc,n) aproxima la solucion en n pasos (en las abscisas x) de:
%        y'' = F.u + F.v y + F.w y' 
%        y(bc.a) = bc.ya
%        y(bc.b) = bc.yb
h = (bc.b - bc.a) / n;
x = bc.a + h *(0:n);
ux = F.u(x); vx = F.v(x); wx = F.w(x);
%matriz:       filas  cols        valor              size 
A =     sparse(1:n-1, 1:n-1,  2 + h^2*vx(2:n),     n-1, n-1); %diag
A = A + sparse(2:n-1, 1:n-2, -1 - h*  wx(3:n)/2,   n-1, n-1); %subdiag
A = A + sparse(1:n-2, 2:n-1, -1 + h*  wx(2:n-1)/2, n-1, n-1); %superdiag
%segundo miembro
b      = -h^2 * ux(2:n);
b(1)   = b(1)   - bc.ya * (-1 - h*wx(2)/2);
b(n-1) = b(n-1) - bc.yb * (-1 + h*wx(n)/2);
%solucion
y = [bc.ya (A\b')' bc.yb];
%2021 Francisco Pena