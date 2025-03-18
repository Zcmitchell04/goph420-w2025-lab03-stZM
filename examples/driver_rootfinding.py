import numpy as np
import matplotlib.pyplot as plt
from src.goph420_lab03.roots import (root_secant_modified,
root_newton_raphson,

                                     )
# Frequency range
frequencies = [0.1, 0.5, 1, 1.4, 2, 5] # Frequencies in Hz

def atote_finding(freq, rho1, rho2, beta1, beta2, H):
    """ Question 2 Part a: uses equation #1

    Parameters:
    ----------

    Returns:
    --------

   """


    nf = len(freq)

    plt.figure()

    zeta_max = np.sqrt(H ** 2 * (beta1 ** -2 - beta2 ** -2))

    for j, f in enumerate(freq):

        def F(z):
            return (rho2/rho1) * np.sqrt(zeta_max**2 - z**2) / z - np.tan(2 * np.pi * f * z)

        atotes = [0.0]
        a = 0.0
        k = 0
        while a < zeta_max:
            a = 0.25 * (2*k+1) / f
            if a < zeta_max:
                atotes.append(a)
            k += 1
        atotes.append(zeta_max)
        n = len(atotes)
        print(n)



        plt.subplot(nf,1,j+1)
        for k, ak in enumerate(atotes):

            if k and k < n-1:
                plt.plot([ak,ak], [-5,5], "--r")

            #plot function F(z)
            if k < n-1:
                zp = np.linspace(ak + 1e-3, atotes[k+1] - 1e-3)

                Fp = F(zp)
                plt.plot(zp, Fp, "-b")

        plt.grid()
        plt.xlabel('zeta')
        plt.ylabel('F(z)')
        plt.title(f'Frequency: {f} Hz')
        plt.ylim(-5,5)
        plt.xlim(0,zeta_max)

    plt.show()





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

    # Calculate wavelengths for each frequency
    lambda_L = c_L / f  # Using equation #3: lambda_L = c_L / f

    # Plot Love wave velocity vs. frequency
    plt.figure(figsize=(12, 6))

    plt.subplot(2, 1, 1)
    plt.plot(f, c_L, label="Love wave velocity (c_L)", color='blue')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Love wave velocity (m/s)')
    plt.title('Love Wave Velocity vs Frequency')
    plt.grid(True)
    plt.legend()

    # Plot Wavelength vs Frequency
    plt.subplot(2, 1, 2)
    plt.plot(f, lambda_L, label="Wavelength (lambda_L)", color='red')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Wavelength (m)')
    plt.title('Wavelength vs Frequency')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


def dispersion_curve(frequencies, rho1, rho2, beta1, beta2, H):
    """ Function to compute Love wave velocity (c_L) and mode solutions for each frequency.
    Parameters:
    frequencies (list): Frequencies to analyze.
    rho1 (float): Density of the surface layer (kg/m^3).
    rho2 (float): Density of the half-space layer (kg/m^3).
    beta1 (float): Shear wave velocity of the surface layer (m/s).
    beta2 (float): Shear wave velocity of the half-space layer (m/s).
    H (float): Thickness of the surface layer (m).
    Returns:
    c_L (list): Love wave velocities for each frequency.
    """
    c_L = []

    # Loop through each frequency to find the corresponding Love wave velocity (c_L)
    for f in frequencies:
        # Define the root-finding problem for each frequency
        def g(zeta):
            return (rho2 / rho1) * (np.sqrt(H ** 2 * (beta1 ** -2 - beta2 ** -2) - zeta ** 2) / zeta) - np.tan(
                2 * np.pi * f * zeta)

        # Solve for zeta using the Newton-Raphson method (or other methods as needed)
        zeta_initial_guess = 0.1  # A reasonable initial guess for zeta
        zeta, _, _ = root_newton_raphson(zeta_initial_guess, g, lambda zeta: (g(zeta + 1e-6) - g(zeta)) / 1e-6)

        # Calculate Love wave velocity c_L for each root of zeta
        c_L_value = H / zeta
        c_L.append(c_L_value)

    return c_L


# Given parameters
rho1 = 1800  # kg/m^3
rho2 = 2500  # kg/m^3
beta1 = 1900  # m/s
beta2 = 3200  # m/s
H = 4000  # m (4 km)


# Optional: Plot asymptotes for root finding for some frequencies
atote_finding(frequencies, rho1, rho2, beta1, beta2, H)

# Find the dispersion curve and plot the results
c_L_values = dispersion_curve(frequencies, rho1, rho2, beta1, beta2, H)

# Plot the mode curves
mode_plt(c_L_values, frequencies)


