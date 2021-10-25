%Enana blanca relativista
hold off
clear all
clc
clf
format long

%Problema de Cauchy
a=1;
B=10^5;
eps0=10^(-3);
c=3*10^8;
alfa=1.473;
beta=52.46;
gamma=4/3;

f=@(r,z) [-alfa*(z(1)^(1/gamma)*z(2))/(r^2),...
          beta*r^2*z(1)^(1/gamma)];

%Pasos
h=10;
n=ceil((B-a)/h);
n=n+1;

m=16;
for i=7:m
%Condición inicial
eta=[10^-i,0];

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