# Utility functions for Lab02

def root_newton_raphson(x0, f, dfdx):
    """This function implements the Newton-Raphson Root Finding Open Method.

    Parameters
    ----------
    x0:     float
            initial guess for value of root
    f:      function
            Function for which root is evaluated
    dfdx:   function
            first derivative of function, f, for which the root is evaluated

    Returns
    -------
    float
    final estimate of the root of the given function

    int
    number of iterations to convergence

    numpy.ndarray
    1D vector of approximate relative error at each iteration
    """