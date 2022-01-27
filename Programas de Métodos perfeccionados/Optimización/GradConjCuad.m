function[xold,nor_grad,index]=GradConjCuad(A,b,xold,itmax,tol,eps)
g=A*xold-b;
d=-g;
alfa=(g.'*g)/(d.'*A*d);
xold=xold+alfa*d;
for i=1:itmax
    g=A*xold-b;
    beta=(d.'*A*g)/(d.'*A*d);
    d=-g+beta*d;
    alfa=(g.'*g)/(d.'*A*d);
    xold=xold+alfa*d;
    
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