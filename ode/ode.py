#RK4
def RK4(t0, x0, tn, n):
    """Método de Runge-Kutta de 4 orden para resolver ecuaciones diferenciales ordinarias.

    Args:
        t0 (float): Primer argumento (el valor de la condición inicial de la variable t)
        x0 (float): Segundo argumento (el valor de la condición inicial de la variable x)
        tn (float): Tercer argumento (el valor de la condición de frontera de la variable t)
        n (int): Cuarto argumento (el número de pasos o iteraciones que se van a utilizar para resolver la ED)

    Examples:
        >>> euler(0, 0, 10, 8)
        (array([ 0.  ,  1.25,  2.5 ,  3.75,  5.  ,  6.25,  7.5 ,  8.75, 10.  ]), array([ 0.        ,  0.60218291,  0.15942164,  0.11425799, -0.77075997,
       -0.16900616,  0.45613153,  0.37341583,  0.23954224]))

    Returns:
        tuple: 2 arreglos, el primero de los n+1 valores de la variable t y el segundo de los n+1 valores de la variable t; en el intervalo [t0,tn]

    """
    h = (tn-t0)/float(n)
    t = np.linspace(t0, tn, n+1)
    x = np.zeros(n+1)
    x[0] = x0
    for i in range(n):
        k1 = h * f(t[i], x[i])
        k2 = h * f(t[i] + 0.5*h, x[i] + 0.5*k1)
        k3 = h * f(t[i] + 0.5*h, x[i] + 0.5*k2)
        k4 = h * f(t[i] + h, x[i] + k3)
        x[i+1] = x[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        x_aprox = x[i+1]
        x_n = x[i]
    return (t, x)


#RK2 
def RK2(t0, x0, tn, n):
    """Método de Runge-Kutta de 2 orden para resolver ecuaciones diferenciales ordinarias.

    Args:
        t0 (float): Primer argumento (el valor de la condición inicial de la variable t)
        x0 (float): Segundo argumento (el valor de la condición inicial de la variable x)
        tn (float): Tercer argumento (el valor de la condición de frontera de la variable t)
        n (int): Cuarto argumento (el número de pasos o iteraciones que se van a utilizar para resolver la ED)

    Examples:
        >>> euler(0, 0, 10, 8)
        (array([ 0.  ,  1.25,  2.5 ,  3.75,  5.  ,  6.25,  7.5 ,  8.75, 10.  ]), array([ 0.        ,  0.73137159,  0.34943373, -0.05274473, -1.14619757,
       -1.26024028, -0.56287143,  0.63843778, -0.11188799]))

    Returns:
        tuple: 2 arreglos, el primero de los n+1 valores de la variable t y el segundo de los n+1 valores de la variable t; en el intervalo [t0,tn]

    
    """
    h = (tn-t0)/float(n)
    t = np.linspace(t0, tn, n+1)
    x = np.zeros(n+1)
    x[0] = x0
    for i in range(n):
        k1 = h * f(t[i], x[i])
        k2 = h * f(t[i] + 0.5*h, x[i] + 0.5*k1)
        x[i+1] = x[i] + k2
        x_aprox = x[i+1]
        x_n = x[i]
    return (t, x)



#Euler 
def euler(t0, x0, tn, n):
    """Método de Euler para resolver ecuaciones diferenciales ordinarias.

    Args:
        t0 (float): Primer argumento (el valor de la condición inicial de la variable t)
        x0 (float): Segundo argumento (el valor de la condición inicial de la variable x)
        tn (float): Tercer argumento (el valor de la condición de frontera de la variable t)
        n (int): Cuarto argumento (el número de pasos o iteraciones que se van a utilizar para resolver la ED)
   
    Examples:
        >>> euler(0, 0, 10, 8)
        (array([ 0.  ,  1.25,  2.5 ,  3.75,  5.  ,  6.25,  7.5 ,  8.75, 10.  ]), array([ 0.        ,  0.        ,  1.18623077, -0.15217513, -0.86222182,
       -1.25962901,  1.19715891,  0.22496448,  0.99163788]))

    Returns:
        tuple: 2 arreglos, el primero de los n+1 valores de la variable t y el segundo de los n+1 valores de la variable t; en el intervalo [t0,tn]
    """
    h = (tn-t0)/float(n)
    t = np.linspace(t0, tn, n+1)
    x = np.zeros(n+1)
    x[0] = x0
    for i in range(n):
        x[i + 1] = x[i] + h*f(t[i], x[i])
    return (t,x)




#Funcion
def f(t, x):
    """La ecuación diferencial (ED) dx/dt a resolver.

    Args:
        t (float): First argument, la variable independiente de la ED
        x (float): Second argument, la variable dependiente de la ED

    Examples:
        >>> f(1.0, 2.0)
        -7.158529015192103

    Returns:
        float: Retorna el resultado de la operación -x**(3) + np.sin(t).

    """
    return (-x**(3) + np.sin(t))







