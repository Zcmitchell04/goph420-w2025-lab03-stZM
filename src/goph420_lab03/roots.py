import numpy as np
import matplotlib

#writitng functions to find the roots using two methods



def root_newton_raphson(x0, f, dfdx):
    """Finding the root using the Newton Raphson method

    Parameters
    --------
    x0 (float) = initial guess
    f (int) = the function
    dfdx (np.ndarray) = the function that gives the first derivative of f

    Returns
    --------
    root_est (float) = gives final estimation of the root
    iterations_needed (int) = The number of iterations needed to converge
    rel_err (np.ndarray) = one dimensional vector, gives approx. relative error AT EACH ITERATION
    """
    x0 =
    f =
    dfdx = np.ndarray()




def root_secant_modified(x0, dx, f) :
    """Finding the root using the Newton Raphson method

        Parameters
        --------
        x0 (float) = initial guess
        dx (int) = the function
        f (np.ndarray) = the function that gives the first derivative of f

        Returns
        --------
        root_est (float) = gives final estimation of the root
        iterations_needed (int) = The number of iterations needed to converge
        rel_err (np.ndarray) = one dimensional vector, gives approx. relative error AT EACH ITERATION
        """
    x0 =
    dx =
    f = np.ndarray()
