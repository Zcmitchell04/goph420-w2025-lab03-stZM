import numpy as np
import matplotlib
from src.goph420_lab03.roots import (root_secant_modified,
root_newton_raphson,

                                     )

# IMPORTANT: if zeta^2 > sqrt(H^2(beta1^-2 - beta2^-2)) then there will be no love waves  (gives zeta_max)
# 2pif(zeta) = (2k+1)(pi/2) --> solve for zeta *f is frequency for this file --> should get (1/4f)(2k+1)
# --> we want these values between 0 and zeta_max
# form root finding problem based on Equation #1 using following assumptions


def root_finding(f, rho1, rho2, beta1, beta2, H):
    """ Question 2 Part a: uses equation #1

    Parameters:
    ----------

    Returns:
    --------

   """

    f = range[2, 50, 0.1]  # change this, just put so it would stop giving me a FREAKING ERROR
    rho1 = 1800  # (kg/m^3)
    rho2 = 2500  # (kg/m^3)
    beta1 = 1900  # (m/s)
    beta2 = 3200  # (m/s)
    H = 4.0  # (km)


def mode_plt(c_L, f):
    """ Question 2 part b: has asymptotes and imports function from Question 1
    -Plots Love wave velocity and zeta(squiggly) over a range of frequencies.
    -Plot a curve for c_L vs frequencies, then each point can determine wavelength lambda_L
    -Plot lambda_L vs. frequencies

    Parameters:
    ----------
    c_L = Love wave velocities
    f = frequencies

    Returns:
    --------
    lambda_L = wavelengths based on c_L vs. Frequencies

    """

    # *If zeta has multiple solutions, the different solutions are called 'modes' and correspond to
    # Love waves of different wavelength lambda_L that are propagating at different velocities c_L; related using Eq #3