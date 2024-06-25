Este sitio contiene la documentación del proyecto *ode-project*, el cual es un módulo para resolver numéricamente ecuaciones diferenciales de primer orden de la forma $\cfrac{dx}{dt} = f(x,t)$ mediante los siguientes tres métodos:

## Método de Euler
Este método se basa en la expansión de Taylor de la función $x(t)$ y se aplica en el punto $t$ para evaluar la derivada y aproximar la función en el punto $x = t + h$.

Para cada punto en el tiempo se aproxima $x$ utilizando el resultado de la iteración anterior avanzando en el tiempo un paso $h$ suficientemente pequeño:

$$\boxed{x(t + h) = x(t) + h \, f(x,t)}$$

El error total de aproximación depende linealmente de $h$ multiplicado por el intervalo en el cual realizamos la integración.
## Método de Runge-Kutta de segundo orden (RK2)
Este método se deriva también de la expansión de Taylor alrededor del punto medio $t + h/2$ aplicando el método de Euler. El paso en cada iteración en este método es $h/2$. 

Las ecuaciones del método RK2 y el punto en el tiempo de la siguiente iteración es:

- $k_1 = h \, f(x,t)$

- $k_2 = h\, f\left(x + \frac{k_1}{2},t + \frac{h}{2}\right)$

$$\boxed{x(t + h) = x(t) + k_2}$$

El error de aproximación de cada paso es de orden $O(h^3)$, mientras que el error global  es de order $O(h^2)$.

## Método de Runge-Kutta de cuarto orden (RK4)
En este método se aplica a más puntos intermedios entre $x(t)$ y $x(t+h)$ por medio de expansiones de Taylor.

Las constantes y el valor de $x$ en la siguiente iteración son:

- $k_1 = h \, f(x, t)$,

- $k_2 = h \, f\left(x + \frac{k_1}{2}, t+\frac{h}2\right)$,

- $k_3 = h \, f\left(x + \frac{k_2}{2}, t+\frac{h}2\right)$,

- $k_4 = h \, f\left(x + k_3, t + h \right)$,

$$\boxed{x(t + h) = x(t) + 1/6 \, (k_1 + 2 k_2 + 2k_3 + k_4)}$$

El error de aproximación es del orden $O(h^5)$, mientras que el error global es aproximadamente del orden $O(h^4)$.
