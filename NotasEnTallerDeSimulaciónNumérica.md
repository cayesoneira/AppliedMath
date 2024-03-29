# Comandos de interés:
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
- Muy interesante la `struct`: una especie de tabla donde pueden introducirse funciones, etc. a las que luego te refieres por etiquetas. Muy potente y útil: se puede usar como argumento de una `function`.
- `nn = [nn norm(nl.yk)]`: básicamente esto es *alargar un vector*.
- Empezamos declarando la función anónima `sol = @(x) exp(x) + cos(x)`, ahora avisamos que *x* es simbólica: `syms x`, acto seguido, sustituyendo en la función anónima la variable simbólica lo que tenemos es precisamente una función simbólica que se puede derivar simbólicamente: `diff(sol(x),x)`, ahora convertimos esta función simbólica en una función anónima que denotamos como *dsol*: `dsol = matlabFunction(diff(sol(x),x))`. Ya podemos crear un vector de ordenadas simplemente haciendo `ye = sol(xh)`.
- Cuando ponemos un `disp` o un `error` lo que ponemos en principio es `error('El algoritmo no converge en iteraciones.')` o quizá `error(nl.maxit)`, es decir, o texto o un número directamente. Si queremos escribir un mensaje que incluya el número y el texto, tenemos que poner **un array de strings**: `error(['El algoritmo no converge en ' num2str(nl.maxit) ' iteraciones.'])`, precisamente con el `[]`. En `title` sí hay que separarlo por comas, en `error` y `disp` parece que no.
- Click en los números de línea en Matlab posibilita crear un punto de stop de correr el código que permite ir ejecutando poco a poco desde ahí e ir viendo qué se ha ido declarando.
- `[a, b, c, d] = deal(1, 2, 3, 4)` sirve para asignar en una sola línea lo que sería `a = 1, b = 2, c = 3, d= 4`.
- `ode45` resuelve un sistema de ecuaciones diferenciales por métodos de Runge Kutta de orden 4 y 5.
- Quizá la forma más cómoda de poner varias figuras juntas: `subplot`.
- `spalloc`hace matrices de ceros, lo que se llaman *matrices sparse* o *matrices dispersas*.
-  `mesher` es una función para crear la malla.
-  `%#ok<NBRAK>`o `%#ok<SPRIX>` evita que salte un aviso de potencial error de Matlab. Tiene interés cuando el autor es consciente de este comportamiento y está seguro de que no es un error. Cambia el nombre según el tipo de aviso.
-  `nnz = 3*m.nnod-2 + (m.deg==2)*2*m.nel;` en esta línea puede verse cómo se usa una condición *if* dentro de la propia línea: si (m.deg==2) es cierto, vale 1, si es falso, vale 0. Y esos números entran en los cálculos. Impresionante esto.
-  `%#ok<SPRIX>` sirve para 
-  `P  = @(x) [1-x, x]; DP = @(x) [ -1, 1];`he aquí la esencia del FEM.
-  `varargin` simplemente refiere a una celda, una *cell*, que es un hueco que se puede llenar con números, letras, char, etc., es parecido a las listas de Python: el hueco más general y que admite todo tipo de clases de elementos, sean números, vectores, strings... Si se usa como argumento de una *function* se espera cualquier cosa: `function u = ef(F, m, bc, qd, varargin)` significa que podemos llamarla como `ef(F, m, bc, qd, @df)` o como `ef(F, m, bc, qd, 2)` o como `ef(F, m, bc, qd, 'patata')`.
-  Las listas etc. se trabajan con llaves `{}`.
-  En la implementación de FEM se puede trabajar con funciones características en la estructura F, pero, a diferencia de hacer subregiones de mallado, esto no permite cambiar el número de mallas en cada región y solo en esa. Esto podría ser interesante si algún material, por ejemplo, requiere por sus características cálculos más precisos.
-  Qué interesante el funcionamiento de los valores lógicos de Matlab: `true` y `false` son valores lógicos y así son interpretados, pero si lo introduces en una operación, por ejemplo `true*2` el resultado es `2`: Matlab lo interpreta como un `1` si lo utilizamos en cálculos.
-  `pdesurf` es capaz, con puntos y aristas como argumento, de saber si le estamos dando una solución por nodos o por vértices al representar soluciones del FEM.
-  Cuando exportamos la solución desde el pdetool obtenemos una "u" que tiene en cada fila el valor de la solución para cada nodo, en orden según la numeración.
-  
---

# Ideas/comentarios:
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
- ¿Por qué ponemos un *or* (||) en lugar de un *and* (&&) en las condiciones de Tolerancia del Residuo y Tolerancia en la Abscisa? Es decir,`if ydif < nl.tol || Fres < nl.tol`. Pues porque si la función es muy muy muy *vaso* entonces converge rápido en la ordenada pero en la abscisa quizá cuando converja el residuo en la ordenada ya sea del orden del épsilon de la máquina.
- El método de Newton para ecs. dif., por contra a iteración funcional, converge muy rápidamente, especialmente con un iterante inicial adecuado. Iteración funcional, sin embargo, perfectamente puede dar lugar a falsas convergencias: converge *a paso de tortuga*.
- **benchmark**: comparativa de rendimiento; básicamente correr un problema resuelto y tabulado con nuestro programa para ver si *tira como debe*.
- Tiene tela que la primera vez en la vida que derivo un operador que incluye operaciones de derivadas (G(y)=y''-u(x)-v(x)y^m-w(x)y') sea en esta asignatura. Bastante peculiar e incomprensible el proceso, de hecho.
- Hay problemas donde la información que se da en la ecuación diferencial es la condición final, no siempre debe ser necesariamente la condición inicial. El hecho es dar tanta información en la variable i-ésima como grados tiene la derivada de dicha variable (i.e. y'' => 2 cond. de contorno).
- La primera vez que tratamos de resolver un problema evolutivo usamos el *método de Euler implícito* para la variable t y *diferencias finitas* para la variable x.
- Para las derivadas primera y segunda hay distintas aproximaciones: podemos hacerlo en i+1 e i, en i-1 e i (de primer orden, lo primero que se ve) o en i-1, i e i+1 (esto son aproximaciones de segundo orden: i-1 e i+1 para la derivada).
- Implementar diferencias finitas hace posible escribir el problema como una ecuación en diferencias (i.e. una ecuación que involucra y_i-1, y_i, y_i+1,, etc.) cuya resolución viene de resolver un sistema matricial (pues al final tenemos un conjunto de ecuaciones lineales).
- Muy interesante este
- Testear código vale dinero. Por eso hay que usarlo lo más posible antes de cambiarlo.
- Ecuaciones diferenciales rígidas (*stiff*) y no rígidas (*non-stiff*): un caso, según Fran, en el que física y métodos numéricos se dan la mano de una forma especial: los *non-stiff* son aquellos que requerirían un paso atento de una iteración a otra por tratarse de modelos de variación muy súbita, y así la física que describen. Los *non-stiff* son mucho más caros por tener que ser más atentos  (por ejemplo con paso variable): los *stiff*, sin embargo, son suficientemente estables y sencillos como para aplicar el mismo método a toda la resolución del problema.
- Aún no alcanzo a comprender por qué, dado que soy neófito en Elementos Finitos, pero al parecer ser puede hacer con polígonos de tantos lados como se quiera.
- El método de diferencias finitas se basa únicamente en plantear un sistema de ecuaciones sobre el vector solución con la información de cómo derivar numéricamente para que converja a la derivada.
- Las condiciones de contorno Robin para el caso de la ecuación del calor tienen perfecto sentido en un contexto de difusión del calor en el borde con dependencia de la temperatura ambiente. Es decir, que, interpretada, esta condición tiene sentido física satisfactorio.

### FEM
- Dice Fran que FEM en varias dimensiones es el mejor recurso en la actualidad para resolver EDPs, de ahí su interés para nuestro estudio.
- COMSOL nació como PDETool de Matlab, dice Fran.
- La gran ventaja de los FEM frente a los método de diferencias finitas es que en 2D y 3D necesitamos mallados tremendamente regulares, pues la ecuación fundamental que permite crear la matriz proviene de una forma de aproximar las derivadas y juega un papel clave el tamaño de los intervalos, etc. En FEM, sin embargo, las regiones son de tamaño arbitrario: esto permite ahorrar una de tiempo de cálculo enorme al centrarnos en las regiones que nos interesan.
- Tratamos de cambiar la ecuación diferencial original por una ecuación integral, lo que se denomina formulación débil, para estar en las hipótesis de Lax-Milgram: un teorema que involucra, de Hilbert, condiciones de existencia de solución para dicho problema. Una especie de *podemos encontrar un elemento ortogonal*. Y esto solo es posible en el lenguaje de integrales, etc., por representar productos escalares o normas.
- The Finite Element Method: Its Basis and Fundamentals by Olek C Zienkiewicz (Author), Robert L Taylor (Author) y J.Z. Zhu (Author). Según Fran, la biblia de los Elementos Finitos.
- **Formas de escribir la condición Dirichlet en FEM.** *Bloqueo por sustitución*: cambia muchísimos elementos de la matriz, pero mantiene la estructura. *Bloqueo por pivote*: cambia el mínimo número de elementos de la matriz, pero puede estropear el número de condicionamiento. *Bloqueo por condensación*: no solo modifica los valores de la matriz, sino también la estructura.
- A partir de FEM cuadrático el número de nodos ya es más grande estrictamente que el número de vértices de la malla.
- Hay dos formas en FEM de *refinar* la resolución: en primer lugar, aumentar el grado; en segundo lugar, estrechar el tamaño del elemento finito, i.e. hacer más fina la malla. Se denominan *métodos hk* a los que van cambiando estos valores de forma óptima.
- Estudiar la matriz del FEM puede ser útil antes de resolver: puede ser útil calcular los autovalores más grandes y quitarle a dicha matriz los cachos que no involucren dichos autovalores. Son la familia de lo *métodos de orden reducido*, muy relacionado también con los métodos de *Singular Value Decomposition (SVD)*.
- Los cálculos numéricos profesionales (ejecución real de FEM en ambiente de empresa, etc.) necesitan un tratamiento local: así es posible repartir el problema entre varios ordenadores y resolverlo así más rápido. Se hace en paralelo el cálculo.
- Naturalmente, no podemos permitirnos perder orden de convergencia del método por una etapa como la integración numérica en FEM. Claro que tampoco podemos permitirnos una calidad inmensa pues la etapa limitante será otra y, al mismo tiempo, no queremos que tarde un tiempo excesivo. Precisamente el interés del tratamiento local del método es, como acabamos de decir, paralelizar y ahorrar tiempo.
- En realidad es bastante impresionante y potente la opción *function* de Matlab: simplemente poder crear tu función de fácil uso y con memoria aislada, etc.
- Fran: *El FEM es el caso final del análisis numérico: todas las disciplinas numéricas modernas surgen a partir de él*.
- En FEM, *x* es la posición global en la malla, la coordenada Euleriana, mientras que *x^* (equis tilde) es la posición relativa dentro del elemento finito, lo que se denomina coordenada lagrangiana.
- Fran: *Hemos llegado al final del juego (el último gran método numérico, FEM), pero no es el final de los juegos*.
- En FEM, la función *F_k* sirve para hacer el cambio de variable desde la x (euleriana) a x tilde (lagrangiana), donde están escritos los polinomios de un elemento finito en general.
- Un pivote, en general, es un valor que se prioriza numéricamente alrededor del cual vamos a calcular. Sea multiplicándolo por un número muy grande para darle prioridad, sea haciendo nulo el resto, sea ignorándolos, etc.
- Hay varios ámbitos en los cuales uno puede optimizar el tiempo: tiempo de programación (usar fsolve, por ejemplo, en lugar de programar un RK4 o algo así) o tiempo de computación (sabiendo que fsolve sale caro computacionalmente, programar un RK4 quizá no sea desmesurado). Todo depende, por supuesto, de qué queramos hacer.
- Fran: *No sabéis cuántos matrimonios ha salvado que se pueda trabajar así de rápido con las matrices sparse en Matlab*.
- La función test *v* de FEM se introduce cuando buscamos escribir el problema en la formulación débil: al integrar por Gauss aparece. Después, con un poco de trajineo, somos capaces de quitarla del medio (al fin y al cabo la formulación débil es para toda función test v).
- Todo el formalismo de la escritura débil de estos problemas viene de la teoría de Espacios de Hilbert.
- En FEM, la formulación débil es para toda función test v: esto significa que también tiene que cumplirse para una cierta v test que elijamos. Aquí es donde entra la aproximación de Galerkin y la forma muy particular de base local polinómica con la que aproximamos dicha v test y que permite escribir el problema de una forma muy concentrada (sparse matrix etc.).
- Para el planteamiento del problema de FEM tenemos dos teoremas clave: Lax-Milgram para existencia y unicidad, Cea para la calidad de la aproximación (el orden, vaya).
- En FEM 2D las esquinas de los dominios, incluso en el cuadrado, no funcionan demasiado bien...
- Normalmente la numeración de los elementos finitos se hace en el sentido de las agujas del reloj.
- Otra fuente de error en FEM: la aproximación del dominio, pues en ocasiones (arcos de circuenferencia, etc.) no es perfecta.
- Tanto triángulos como rectángulos se pueden usar para mallar.
- Dos elementos finitos pueden tener en común: caras (solo 3D), aristas (solo 2D), vértices o nada. No es posible que estén a caballo entre estas opciones. Y esto es una triangulación como hacemos en Topología de Superficies tal cual.
- Mallas adaptativas: cuanto más grande es el gradiente en una región, más pequeño interesa hacer el elemento finito para que la norma de Sobolev siga siendo pequeña.
- Dicen los apuntes: **Si una geometría tiene varios materiales se recomienda que cada elemento quede completamente contenido en uno de ellos**.
- Hay formas de medir la calidad de la malla, como por ejemplo, dentro de un elemento finito, la razón entre radios de las circunferencias inscritas y circunscritas. Un valor muy alto de ese cociente señala la cercanía de un elemento finito a ser degenerado.
- Fases del trabajo de modelado: preproceso (mallado, elección del sistema y condiciones de contorno, etc.), simulación (resolución), posproceso (análisis y visualización de resultados).
- Más importante que replicar los resultados de un cierto benchmark (que dependen en algunos detalles del propio PC), lo importante a la hora de ver la implementación de un algoritmo es que el orden del mismo se corresponda, es decir, que se replique la gráfica de orden vs anchura de los elementos.
- El épsilon de la máquina puede compararse literalmente como un ruido blanco, un fondo donde no se puede trabajar exactamente con ciertos valores.
- De nuevo, un caso donde el límite de la matemática se ve interrumpido por el límite numérico.
- Al parecer, es un convenio a nivel mundial escribir un problema de valor inicial escribiendo la expresión de la ecuación *en* (*in*) el dominio pero la expresión de las condiciones de contorno *sobre* (*on*) la frontera. Peregrina lo dijo con bastante seriedad.
- En el método de Newton, cuando queremos tomar un valor para la siguiente iteración, interpolamos (linealmente en FEM de grado 1 y con parábolas en FEM de grado) en lugar de tomar el valor más cercano o métodos más complejos porque así cometemos el mismo error que el propio FEM.
- Cuando aproximamos  en H1(0,1), lo mejor que conseguimos en general es O(h^2).
- Si la ecuación diferencial es simétrica, la matriz elemental es simétrica.
- O(h^2), es decir, *'O' grande*, significa que el límite del cociente con h^2 es un valor real, pero no necesariamente nulo. Eso es o(h^2), es decir, *'o' pequeña*.
