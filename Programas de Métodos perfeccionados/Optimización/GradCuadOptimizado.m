function[xold,nor_grad,index]=GradCuadOptimizado(A,b,xold,itmax,tol)
gradJ=@(x) A*x-b;
g=gradJ(xold);
for i=1:itmax
    
    aux=A*g;
    alfa=-(g.'*g)/(g.'*aux);
    xold=xold+alfa*g;
    
    if norm(g,2)<tol
        index=i;
        nor_grad=norm(g,2);
        break;
    else
        index=-1;
        nor_grad=Inf;
    end
    g=g+alfa*aux;
    
end