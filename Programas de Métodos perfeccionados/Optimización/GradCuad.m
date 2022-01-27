function[xold,nor_grad,index]=GradCuad(A,b,xold,itmax,tol)
gradJ=@(x) A*x-b;
for i=1:itmax
    alfa=-(gradJ(xold).'*gradJ(xold))/(gradJ(xold).'*A*gradJ(xold));
    xold=xold+alfa*gradJ(xold);
    if norm(gradJ(xold),2)<tol
        index=i;
        nor_grad=norm(gradJ(xold),2);
        break;
    else
        index=-1;
        nor_grad=Inf;
    end
end
end

