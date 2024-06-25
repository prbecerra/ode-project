# Numerical methods to solve one-variable ODEs
def euler(x0,t,f):
    """ Compute and return the solution of a first-order one-variable differential equation (ODE) by Euler's method.

    Examples:
        Solving dx/dt = -x^3 + sin(t)
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> def f(x,t): return -x**3 + np.sin(t)
        >>> x0 = 0
        >>> t = np.linspace(0,10,100)
        >>> x_euler(x0,t,f)
        >>> plt.plot(t,x_euler)
        >>> plt.show
    
    Args:
        x0 (float): First argument. A floating number representing the initial condition.
        t (array): Second argument. A one-dimensional array of time with the start and end of the interval and the number of steps.
        f (function): Third argument. The function to solve which form is: dx/dt = f(x,t)

    Returns:
        x(array): A one-dimensional array with the aproximate solution of dx/dt = f(x,t).
    """
    h = t[1] - t[0]                 # Size of a single step
    x = np.zeros(t.size)            # One-dimensional array initialized to zero of size t
    x[0] = x0                       # Initial condition
    
    for i in range(t.size-1):
        x[i+1] = x[i] + (h * f(x[i],t[i]))
    return x

def rk2(x0,t,f):
    """Compute and return the solution of a first-order one-variable differential equation (ODE) by Ruge-Kutta second order's method.

    Examples:
        Solving dx/dt = -x^3 + sin(t)
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> def f(x,t): return -x**3 + np.sin(t)
        >>> x0 = 0
        >>> t = np.linspace(0,10,100)
        >>> x_rk2(x0,t,f)
        >>> plt.plot(t,x_rk2)
        >>> plt.show
    
    Args:
        x0 (float): First argument. A floating number representing the initial condition.
        t (array): Second argument. A one-dimensional array of time with the start and end of the interval and the number of steps.
        f (function): Third argument. The function to solve which form is: dx/dt = f(x,t)

    Returns:
        x (array): A one-dimensional array with the aproximate solution of dx/dt = f(x,t).
    """
    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i] + k1/2 , t[i] + k1/2)
        x[i+1] = x[i] + k2
    return x

def rk4(x0,t,f):
    """Compute and return the solution of a first-order one-variable differential equation (ODE) by Ruge-Kutta fourth order's method.

    Examples:
        Solving dx/dt = -x^3 + sin(t)
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> def f(x,t): return -x**3 + np.sin(t)
        >>> x0 = 0
        >>> t = np.linspace(0,10,100)
        >>> x_rk4(x0,t,f)
        >>> plt.plot(t,x_rk4)
        >>> plt.show
    
    Args:
        x0 (float): First argument. A floating number representing the initial condition.
        t (array): Second argument. A one-dimensional array of time with the start and end of the interval and the number of steps.
        f (function): Third argument. The function to solve which form is: dx/dt = f(x,t)

    Returns:
        x (array): A one-dimensional array with the aproximate solution of dx/dt = f(x,t).
    """
    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i] + k1/2 , t[i] + k1/2)
        k3 = h * f(x[i] + k2/2 , t[i] + k2/2)
        k4 = h * f(x[i] + k3 , t[i] + h)
        x[i+1] = x[i] + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x

