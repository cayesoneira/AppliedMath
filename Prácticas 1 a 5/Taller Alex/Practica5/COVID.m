% PRÁCTICA 5 COVID

% lo que está en código nos permite obtener los resultados del paper

clear all
global tau beta delta gamma eta theta mu nu sigma Lambda r_1 r_2

% Parámetros del modelo de refencia [1]
tau=0.0002; beta=0.0805; delta=1.6728e-5;
gamma=2.0138e-4; eta=0.4478;
theta=0.0101; mu= 0.0106;
nu=3.2084e-4; sigma=0.0668;
Lambda=0.02537; r_1=5.7341e-5; r_2=1.6728e-5;


%  DESCOMENTAR PARA CONSIDERAR UN NUEVO MODELO MÁS COHERENTE
%     % Parámetros mdificados del modelo de refencia [1]
%     % El objetivo ha sido considerar una tasa de natalidad y de mortalidad
%     % y unos periodos de recuperación razonables
%
%     Lambda=Lambda/365; beta=beta*365; mu=mu/365;
%     r_1=5.7341e-5*1000; r_2=1.6728e-5*1000;
% 
%     %Incremento del número básico de reproducción
%     beta=beta*100;


% Situación inicial
S_0=0.5; E_0=0.2; Q_0=0.1; I_A_0=0.2; I_S_0=0.1; R_0=0;
N_0=S_0+E_0+Q_0+I_A_0+I_S_0;
disp(['La población inicial es ', num2str(N_0),' individuos'])

% CÁLCULO DEL NÚMERO BÁSICO DE REPRODUCCIÓN
R0=beta*Lambda/((gamma+mu+eta+sigma)*(tau+mu));
disp(['El número básico de reproducción para los datos Tabla es ', num2str(R0),' días^(-1)'])

% Si el número de individuos es Pobl
% Pobl=1000;    
% k=Pobl/N_0;
% Lambda_k=Lambda*k;
% beta_k=beta/k;
% R0_k=beta_k*Lambda_k/((gamma+mu+eta+sigma)*(tau+mu));
% disp(['El número básico de reproducción para ', num2str(Pobl),' individuos es ', num2str(R0_k),' días^(-1)'])


% CÁLCULO DEL NÚMERO DE DÍAS PREVISTO PARA LA RECUPERACIÓN PROMEDIO DE UN ASINTOMÁTICO
RI_A=1/r_1;
disp(['El número de días previsto para la recuperación de un asintomático es ', num2str(RI_A), ' días'])


% CÁLCULO DEL NÚMERO DE DÍAS PREVISTO PARA LA RECUPERACIÓN PROMEDIO DE UN SINTOMÁTICO
RI_S=1/r_2;
disp(['El número de días previsto para la recuperación de un sintomático es ', num2str(RI_S), ' días'])


% CÁLCULO DEL NÚMERO DE DÍAS PREVISTO PARA LA INCUBACIÓN DE LA ENFERMEDAD
DInc=1/(nu+sigma);
disp(['El número de días previstos de incubación es ', num2str(DInc), ' días'])


% CÁLCULO DE LA TASA DE NATALIDAD POR CADA MIL HABITANTES
TNat=Lambda*365*1000;
disp(['Tasa de natalidad al año por mil habitantes ', num2str(TNat),'% al año'])
disp(['Tasa de natalidad al año por cien habitantes ', num2str(TNat/10),'% al año'])


% CÁLCULO DE LA TASA DE MORTALIDAD POR CADA MIL HABITANTES
TMort=mu*365*1000;
disp(['Tasa de mortalidad al año por mil habitantes ', num2str(TMort),'% al año'])
disp(['Tasa de mortalidad al año por cien habitantes ', num2str(TMort/10),'% al año'])

% FIGURA 5.5 DE [1]
opciones=odeset('RelTol',1e-7);
[t,y]=ode45('f',[0,100],[0.5, 0.2, 0.1, 0.2, 0.1, 0]', opciones);
c='grbcy*';
titulo=['Susceptibles Expuestos Cuarentena Asintomáticos Sintomáticos Retirados'];
l=[0 12 22 33 47 60 70];
figure(1)
hold on 
for i=1:2:5
    subplot(3,2,i)
        plot(t,y(:,i),c(i))
        title(titulo(l(i)+1:l(i+1)))
    xlabel('Semanas')
    subplot(3,2,i+1)
        plot(t,y(:,i+1),c(i+1))
        title(titulo(l(i+1)+1:l(i+2)))
        xlabel('Semanas')
end
hold off

% Gráfica corregida considerando 100 días en vez de 100 semanas
figure(2)
hold on 
for i=1:2:5
    subplot(3,2,i)
        plot(t,y(:,i),c(i))
        title(titulo(l(i)+1:l(i+1)))
    xlabel('Días')
    subplot(3,2,i+1)
        plot(t,y(:,i+1),c(i+1))
        title(titulo(l(i+1)+1:l(i+2)))
        xlabel('Días')
end
hold off     




