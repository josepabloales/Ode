#Tarea 4: Documentación de código y Git
## Estudiantes: José Pablo Alvarado, Mateo Heinz, Ana Laura Víquez


### Métodos de resolución numérica de ecuaciones diferenciales ordinarias
### 1) Método de Euler

El método de Euler se basa en la expansión de Taylor de la función $x(t)$:
$$
\text{Expansión Taylor} \Rightarrow x(t+h) = x(t) + h\frac{dx}{dt} + \overbrace{ \frac{h^2}{2} \frac{d^2x}{dt^2} } ^{\epsilon} + O(h^3).
$$
Se evalua la función x en un tiempo t más un $h$ suficientemente pequeño:
$$
\boxed{x(t + h) = x(t) + hf(x,t).}
$$
Esta ecuación da una aproximación relacionada con el número de pasos en el tiempo que utilicemos para llegar a ella y el error asociado con la aproximación también esta determinado por este número de pasos como se muestra en la siguiente ecuación:

$$
\sum\epsilon = \sum_{k=0}^{N-1}\frac{h^2}{2}\left. \frac{d^2x}{dt^2} \right|_{x_k, t_k} = \frac{h}{2}\sum_{k=0}^{N-1}h\left.\frac{df}{dt}\right|_{x_k, t_k}\\
$$

$$
\approx \frac{h}2\int_a^b\frac{df}{dt}d t = \frac{h}{2}\left[f_b - f_a\right].
$$
En la ecuación anterior se asume que se toma $N = (b-a)/h$ pasos temporales para llegar al punto final.

El error total de aproximación depende $h$ linealmente multiplicado por el intervalo en el cual realizamos la integración.

* Para algunas aplicaciones, esto es suficiente. Para otras, necesitamos una mejor aproximación.
* El algoritmo toma la siguiente forma:
  - Empezar con $t = t_0$, $x = x_0$
  - Discretizar el tiempo en pasos temporales de forma equidistante con espaciamiento $h$, donde cada punto en el tiempo está denotado con $t_i$
  - Para cada punto en el tiempo encontrar $x$ utilizando el resultado de la iteración previa: $x_i = x_{i-1} + hf(x_{i-1})$





### 2) Método de Runge-Kutta


#### Rk2
La idea del método RK2 es utilizar el punto medio para evaluar el método de Euler, como se indica en la figura. Mientras que el método de Euler se aplica en el punto $t$ para evaluar la derivada para aproximar la función en el punto $x = t + h$, el método RK2 utiliza el punto medio $t + h/2$. 

De esta forma, se alcanza una mejor aproximación para el mismo valor de $h$.

El método se deriva aplicando la serie de Taylor alrededor del punto medio $t + h/2$ para obtener el valor de la función en el punto $x(t + h)$. Tenemos:

$$
x(t + h) = x\left(t + \frac{h}{2}\right) + \frac{h}{2}\left(\frac{{\rm d}x}{{\rm d}t}\right)_{t+h/2} + \frac{h^2}{8}\left(\frac{{\rm d}^2x}{{\rm d}t^2}\right)_{t+h/2} + O(h^3)
$$


Similarmente, podemos hacer lo mismo para $x(t)$, tal que:

$$
x(t) = x\left(t + \frac{h}{2}\right) - \frac{h}{2}\left(\frac{{\rm d}x}{{\rm d}t}\right)_{t+h/2} + \frac{h^2}{8}\left(\frac{{\rm d}^2x}{{\rm d}t^2}\right)_{t+h/2} + O(h^3)
$$


Al sustraer ambas ecuaciones obtenemos
$$
x(t + h) = x(t) + h\left(\frac{{\rm d}x}{{\rm d}t}\right)_{t+h/2} + O(h^3)
$$

Finalmente,
$$
\boxed{x(t + h) = x(t) + hf[x(t + h/2), t + h/2] + O(h^3)}.
$$

El término de orden $h^2$ desaparece y nuestra aproximación tiene un error de orden $h^3$. Recordemos que incrementar el orden del error por un orden de magnitud es muy beneficioso a nivel computacional. 

El único problema es que requerimos conocer el valor de la función en el punto medio $x(t + h/2)$, el cual desconocemos.

Para aproximar este valor utilizamos el método de Euler con un paso $h/2$, $(x + h/2) = x(t) + \frac{h}{2}f(x,t)$. De esta manera, obtenemos las ecuaciones del método RK2:

* $k_1 = hf(x,t),$
* $k_2 = hf\left(x + \frac{k_1}{2},t + \frac{h}{2}\right)$
* $x(t + h) = x(t) + k_2$
El error de aproximación de cada paso es de orden $O(h^3)$, mientras que el error global (con un análisis similar al que hicimos con el método de Euler) es de order $O(h^2)$. 

Cabe recalcar que al utilizar el método de Euler para la primera parte de la aproximación, el error también es de $O(h^3)$ y por ende el error de aproximación se mantiene de $O(h^3)$.






#### Rk4
Este método es igual al RK2 pero se toman en vez de 2 valores $h_1$ y $h_2$, 4 valores $h_1$, $h_2$, $h_3$ y $h_4$. Y estos valores determinan la ecuación $x(t+h)$ de una mejor manera, haciendo la aproximación más exacta y convergiendo más rápido al resultado real de la ecuación diferencial.

Los valores de $h_n$ son:

* $k_1 = hf(x, t)$,
* $k_2 = hf\left(x + \frac{k_1}{2}, t+\frac{h}2\right)$,
* $k_3 = hf\left(x + \frac{k_2}{2}, t+\frac{h}2\right)$,
* $k_4 = hf\left(x + k_3, t + h \right)$,
* $x(t+h) = x(t) + \frac{1}{6}(k_1 + 2 k_2 + 2k_3 + k_4)$.

Generalmente este es el mejor método para programar soluciones a ecuaciones diferenciales por su facilidad de programación y la precision de sus resultados. 

El error de aproximación es $O(h^5)$, mientras que el error global es aproximadamente del orden $O(h^4)$.