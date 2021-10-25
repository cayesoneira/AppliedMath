function[xold,nor_grad,index,normas]=GradConjCuadOptimizado(A,b,xold,itmax,tol,eps)
g=A*xold-b;
d=-g;
aux=A*d;
alfa=(g.'*g)/(d.'*aux);
xold=xold+alfa*d;
for i=1:itmax
    beta=((alfa*aux).'*(g+alfa*aux))/(g.'*g);
    g=g+alfa*aux;
    d=-g+beta*d;
    aux=A*d;
    alfa=(g.'*g)/(d.'*aux);
    xold=xold+alfa*d;
    normas(i)=norm(g,2);
    
    if norm(g,2)<tol
        index=i;
        nor_grad=norm(g,2);
        break;
    else
        index=-1;
        nor_grad=Inf;
    end
    
    if beta<eps
        REVENTASTELAMAQUINAWEY=1
        break;
    end
    if alfa<eps
        REVENTASTELAMAQUINAWEY=1
        break;
    end
end