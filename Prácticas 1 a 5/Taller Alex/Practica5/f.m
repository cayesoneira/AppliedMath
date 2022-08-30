function z=f(t,y)
% modelo SEQIR evolución COVID-19
z=zeros(6,1);

global tau beta delta gamma eta theta mu nu sigma Lambda r_1 r_2
z(1)=Lambda-(tau+mu)*y(1)-beta*y(1)*y(2);
z(2)=beta*y(1)*y(2)-(gamma+mu+eta+sigma)*y(2);
z(3)=tau*y(1)+gamma*y(2)-(mu+nu+theta)*y(3);
z(4)=sigma*y(2)+theta*y(3)-(mu+r_1)*y(4);
z(5)=eta*y(2)+nu*y(3)-(delta+mu+r_2)*y(5);
z(6)=r_1*y(4)+r_2*y(5)-mu*y(6);