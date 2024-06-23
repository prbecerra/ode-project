
def euler(x0,t,f):
    """ Compute and return the solution of a first-order one-variable differential equation by Euler's method.

    Examples:
        >>> euler(x0,t,f)
        x

    Args:
        x0 (float): First argument. A floating number representing the initial condition.
        t (vector): Second argument. A one-dimensional array of time.
        f : Third argument. The function to solve which form is: dx/dt = f(x,t)

    Returns:
        x: A one-dimensional array with the solution of dx/dt = f(x,t).
    """
    
    h = t[1] - t[0]                 # Size of a single step
    x = np.zeros(t.size)            # One-dimensional array initialized to zero of size t
    x[0] = x0                       # Initial condition
    
    for i in range(t.size-1):
        x[i+1] = x[i] + (h * f(x[i],t[i]))
    return x


def rk2(x0,t,f):
    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i] + k1/2 , t[i] + k1/2)
        x[i+1] = x[i] + k2
    return x

def rk4(x0,t,f):
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


def rk4(func, ham, y_n, h):
    k1 = func(ham,y_n)
    k2 = func(ham,y_n+h/2*k1)
    k3 = func(ham,y_n+h/2*k2)
    k4 = func(ham,y_n+h*k3)
    
    return y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    # Esta función debe devolver y_{n+1}

    
