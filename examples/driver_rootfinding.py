import numpy as np
import matplotlib as plt
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

    f = [0.1, 0.5, 1.0, 2.0]  # change this, just put so it would stop giving me a FREAKING ERROR
    nf = len(f)

    plt.figure()

    rho1 = 1800  # (kg/m^3)
    rho2 = 2500  # (kg/m^3)
    beta1 = 1900  # (m/s)
    beta2 = 3200  # (m/s)
    H = 4.0  # (km)

    asymptote = np.sqrt(H ** 2 * (beta1 ** -2 - beta2 ** -2))
    atotes = [0.0]
    a = 0.0
    k = 0

    zeta_max = H ** 2 * (beta1 ** -2 - beta2 ** -2)

    for j, f in enumerate(f):

        def F(z):
            return (rho2/rho1) * np.sqrt(zeta_max**2 - z**2) / z - np.tan(2 * np.pi * f * z)
        while a < zeta_max:
            a = 0.25 * (2*k+1) / f
            if a < zeta_max:
                atotes.append(a)
            k += 1
        atotes.append(zeta_max)
        n = len(atotes)


        plt.subplot(nf,1,j+1)
        for k, ak in enumerate(atotes):

            if k and k < n-1:
                plt.plot([ak,ak], [-5,5], "--r")

            if k < n-1:
                zp = np.linspace(ak + 1e-3, atotes[k+1] - 1e-3)

                Fp = F(zp)
                plt.plot(zp, Fp, "--b")

        plt.grid()



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