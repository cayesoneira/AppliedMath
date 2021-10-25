%Estrella de neutrones pura no relativista
hold off
clear all
clc
clf
format long

%Problema de Cauchy
a=1;
B=200000;
eps0=0.08969;
c=3*10^8;
alfa=1;
beta=0.7636;
K=1.914;
gamma=5/3;

f=@(r,z) [-1/(r^2)*alfa*z(2)^(1/gamma)*(1+K^(1/gamma)*z(1)^(1-1/gamma))*...
    (1+4*pi*eps0*(r^3*z(1))/(z(2)))*(1-2.94*z(2)/r)^(-1),...
    beta*r^2*z(1)^(1/gamma)];

%Pasos
h=10;
n=ceil((B-a)/h);
n=n+1;

m=10;
for i=5:m
%Condición inicial
eta=[1,1];

%Funciones de los métodos particulares
tic
 [x,y,p]=classicrkfunc(f,a,eta,h,n);
% [x,y,p]=gausslegendrerk2func(f,a,eta,h,n);
toc

rlist(i)=x(p-1);
mlist(i)=y(p-1,2);

hold on
subplot(1,3,1),plot(log10(x),y(:,1))
subplot(1,3,2),plot(log10(x),y(:,2))
end

subplot(1,3,3),plot(log10(rlist),mlist)