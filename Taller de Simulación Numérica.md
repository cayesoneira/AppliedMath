### Comandos de interés:

- `errorbar(abscisa,ordenada,vector de errores);`
- `./` divide componente a componente: es el efecto del punto delante del signo de operación
- `plot(h(1,end-1),p)` si p mide una posición más que h; cogemos hasta la penúltima componente de h
- `axis` va con cuatro valores: abscisa mínima, abscisa máxima, ordenada mínima, ordenada máxima.
- la derivada simbólica en matlab no necesita especificar la variable si solo hay una variable; de hecho, matlab da prioridad a la x.
- `contour` hace las curvas de nivel; `contourf` las *rellena* (fill)
- `repmat`: *repite un vector para hacer una matriz*
- `trapz` es aplicar la regla del trapecio dando un vector de abscisas y otro de ordenadas.
- `orden(F, bc, sol, @normaInf, 'Orden de convergencia en norma L^\infty', 1.9, 2.1)`: esto es un ejemplo de mi *function* que calcula el orden. Lo llamativo es que uno de los argumentos de mi *function* es otra *function*, lo que hay que hacer es poner un `@` para convertirla en una **función anónima**. Si el argumento es directamente una función anónima (por ejemplo declarada tipo `f=@(x) x^2`) entonces no hace falta poner el `@` cuando se usa como argumento.
- Formatos de número en `disp`: `disp([' Dx = ' num2str(dx,'%4.3f')])` dice con el `%4.3f` que quiere un *float* con ciertos decimales, etc.
- 
---
### Ideas/comentarios:

- Condiciones de contorno Dirichlet (fijar la solución), Neumann (fijar las derivadas de la solución), Robin (combinación lineal de los valores de una función y los valores de su derivada sobre la frontera del dominio), mixtas[?] (combinación cualquiera)
- Por ahora la asignatura parece ir poco a poco, explicando matlab desde la base (las clases son programación guiada con la explicación de los comandos).
- para un pc de 64 bits el eps de la maq es `2.22e-16`
- Matlab se lleva mal con la memoria dinámica: es mejor que los vectores se declaren directamente con su dimensión. Lo otro podría causar que la memoria deba estar realojándose, lo cual ralentiza el proceso.
- la matemática computacional tiene dos límites: el matemático (h's pequeños deben usarse) y el computacional (h's más grandes que el eps de la máquina); hay un intervalo para el h con el que debemos trabajar.
- Concepto de *test académico*: eufemismo elegante para decir *ejercicio*
- Método de elementos finitos >> Método de diferencias finitas >> Métodos para problemas de valor inicial >> Runge Kutta's
- Tener condiciones Dirichlet se denomina también *bloquear*; cond dirichlet en cierto extremo se dice *bloquear* dicho extremo.
- Con condiciones Dirichlet los errores crecen cuando nos acercamos al centro del intervalo: en los extremos ya estamos dando la solución. No pasa lo mismo con condiciones mixtas, que pueden ser fuente de error de por sí.
- Evolutivo = No estacionario.
- El numérico aporta intuición sobre qué temas teóricos pueden tener solución o no (temas de existencia y unicidad, etc.) y si perder el tiempo en ellos.
- La iteración funcional se puede aplicar también a EDs considerando un operador que incluya las acciones de *derivar*, entre otras operaciones. El objetivo de la iteración funcional siempre ha sido escribir la solución como un punto fijo de una aplicación. El tema 4 hace un interesante ejercicio sobre estas ideas.
- En el núcleo del algoritmo originario de Deep Learning se encuentra un algoritmo de Newton, esto es, una aproximación de un problema por otro problema lineal.
- ¿Por qué ponemos un *o* (||) en lugar de un *y* (&&) en las condiciones de Tolerancia del Residuo y Tolerancia en la Abscisa?
> `if ydif < nl.tol || Fres < nl.tol`
Pues porque si la función es muy muy muy *vaso* entonces converge rápido en la ordenada pero en la abscisa quizá cuando converja el residuo en la ordenada ya sea del orden del épsilon de la máquina.
