### Comandos de interés:

- `errorbar(abscisa,ordenada,vector de errores);`
- `./` divide componente a componente: es el efecto del punto delante del signo de operación
- `plot(h(1,end-1),p)` si p mide una posición más que h; cogemos hasta la penúltima componente de h
- `axis` va con cuatro valores: abscisa mínima, abscisa máxima, ordenada mínima, ordenada máxima.
- la derivada simbólica en matlab no necesita especificar la variable si solo hay una variable; de hecho, matlab da prioridad a la x.
- `contour` hace las curvas de nivel; `contourf` las *rellena* (fill)
- `repmat`: *repite un vector para hacer una matriz*
- `trapz` es aplicar la regla del trapecio dando un vector de abscisas y otro de ordenadas.
---
### Ideas/comentarios:

- Condiciones de contorno Dirichlet (fijar la solución), Neumann (fijar las derivadas de la solución), Robin (combinación lineal de los valores de una función y los valores de su derivada sobre la frontera del dominio), mixtas[?] (combinación cualquiera)
- Por ahora la asignatura parece ir poco a poco, explicando matlab desde la base (las clases son programación guiada con la explicación de los comandos).
- para un pc de 64 bits el eps de la maq es `2.22e-16`
- Matlab se lleva mal con la memoria dinámica: es mejor que los vectores se declaren directamente con su dimensión. Lo otro podría causar que la memoria deba estar realojándose, lo cual ralentiza el proceso.
- la matemática computacional tiene dos límites: el matemático (h's pequeños deben usarse) y el computacional (h's más grandes que el eps de la máquina); hay un intervalo para el h con el que debemos trabajar.
- Concepto de *test académico*: eufemismo elegante para decir *ejercicio*
- Método de elementos finitos >> Método de diferencias finitas >> Métodos para problemas de valor inicial >> Runge Kutta's
- Tener condiciones Dirichlet se denomina también *bloquear*; cond dirichlet en cierto extremo 
