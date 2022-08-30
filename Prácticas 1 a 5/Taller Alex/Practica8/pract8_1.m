% PRÁCTICA 8 ejercicio 1
% Ejecutar primero la función Practica8_1.m
% Exporta a la ventana de comandos la MALLA y la SOLUCIÓN

% cálculo del máximo error absoluto
x=p(1,:)'; y=p(2,:)'; norm(u+3*(x.^2+y.^2)/2,Inf)