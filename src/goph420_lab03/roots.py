import numpy as np
import matplotlib

#writitng functions to find the roots using two methods



def root_newton_raphson(x0, f, dfdx, tol = 1e-8, max_iter = 100):
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
    x = x0
    errors = []

    for iter_num in range(1, max_iter + 1):
        # Compute the next approximation using Newton-Raphson formula
        fx = f(x)
        dfx = dfdx(x)

        if dfx == 0:  # Avoid division by zero
            raise ValueError("Derivative is zero, unable to proceed with Newton-Raphson.")

        x_new = x - fx / dfx
        error = abs(x_new - x) / abs(x_new)
        errors.append(error)

        # Check for convergence
        if error < tol:
            return x_new, iter_num, np.array(errors)

        x = x_new

    raise ValueError("Maximum number of iterations exceeded.")



def root_secant_modified(x0, dx, f, tol = 1e-8, max_iter = 100) :
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

    x = x0
    errors = []

    for iter_num in range(1, max_iter + 1):
        # Approximate the derivative using the Modified Secant method
        fx = f(x)
        fx_dx = f(x + dx)

        dfx = (fx_dx - fx) / dx
        if dfx == 0:  # Avoid division by zero
            raise ValueError("Derivative approximation is zero, unable to proceed.")

        x_new = x - fx / dfx
        error = abs(x_new - x) / abs(x_new)
        errors.append(error)

        # Check for convergence
        if error < tol:
            return x_new, iter_num, np.array(errors)

        x = x_new

    raise ValueError("Maximum number of iterations exceeded.")
