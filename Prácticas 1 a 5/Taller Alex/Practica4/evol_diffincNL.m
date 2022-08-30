function [x,t,Y]=evol_diffincNL(l,u,v,w,y0,a,b,npas,t0,T,npast,datCC,datNL)

Dt = (T-t0)/npast;  % incremento de tiempo
t = t0+[0:npast]*Dt;    
h = (b-a)/npas;     % incremento de espacio
x = a+[0:npas]*h;
Y = sparse(npas+1,npast+1); % matriz solucion (npas+1) x (npast+1); toda ceros, la primera columna es y0
Y(:,1)=y0(x);

if datNL.p==1
    % caso LINEAL
    for i=2:length(t)
        u1= @(x) u(x,t(i),datNL.p)-l(x,t(i)).*Y(:,i-1)'/Dt;
        v1= @(x) v(x,t(i))+l(x,t(i))/Dt;
        w1= @(x) w(x,t(i));
        alfa1= datCC.alfa(t(i));
        beta1= datCC.beta(t(i));
        [x,Y(:,i)] = diffinc(u1,v1,w1,a,b,alfa1,beta1,npas,datCC.iopc,datCC.gamma,datCC.delta);
    end
    
% OJO!! HAY QUE REDEFINIR LAS FUNCIONES u,v,w,l PORQUE CAMBIA LA SOLUCION
% DEL PROBLEMA PARA p!=1
elseif datNL.iopNL=='IF' 
    % caso ITERACIÓN FUNCIONAL
     for j=2:length(t)
        up_IF= @(x) u(x,t(j),datNL.p)-l(x,t(j)).*Y(:,j-1)'/Dt;
        vp_IF= @(x) v(x,t(j));
        wp_IF= @(x) w(x,t(j));
        lp_IF= @(x) l(x,t(j))/Dt;
        alfap_IF= datCC.alfa(t(j));
        betap_IF= datCC.beta(t(j));
        [x,Y(:,j),k] = diffincNoLinealIF_evol(lp_IF,up_IF,vp_IF,wp_IF,datNL.p,a,b,alfap_IF,betap_IF,...
            npas,datCC.iopc,datCC.gamma,datCC.delta,datNL.tol,datNL.maxit);
        if k==datNL.maxit
            disp(['No hay convergencia para el tiempo t=',num2str(t0+j*Dt)])
        end
     end
     
elseif datNL.iopNL=='NW' 
    % caso NEWTON
    for j=2:length(t)
        up_N= @(x) u(x,t(j),datNL.p)-l(x,t(j)).*Y(:,j-1)'/Dt;
        vp_N= @(x) v(x,t(j));
        wp_N= @(x) w(x,t(j));
        lp_N= @(x) l(x,t(j))/Dt;
        alfap_N= datCC.alfa(t(j));
        betap_N= datCC.beta(t(j));
        [x,Y(:,j),k] = diffincNoLinealN_evol(lp_N,up_N,vp_N,wp_N,datNL.p,a,b,alfap_N,betap_N,...
            npas,datCC.iopc,datCC.gamma,datCC.delta,datNL.tol,datNL.maxit);
        if k==datNL.maxit
            disp(['No hay convergencia para el tiempo t=',num2str(t0+j*Dt)])
        end
    end
    
end
