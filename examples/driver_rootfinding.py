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

    plt.savefig("C:/Users/zcmit/git/Goph420projects/goph420-w2025-lab03-stZM/figures/atote-plot.png", dpi=300)
    plt.show()



def dispersion_curve(frequencies, rho1, rho2, beta1, beta2, H):
    """ Function to compute Love wave velocity (c_L) and mode solutions for each frequency. """

    mode = 'newton'
    c_L_all_modes = []  # Store Love wave velocities for each mode at each frequency
    zeta_all_modes = []  # Store all the zeta values (modes) for each frequency

    for f in frequencies:
        # Define the root-finding problem for each frequency
        def g(zeta):
            return (rho2 / rho1) * (np.sqrt(H ** 2 * (beta1 ** -2 - beta2 ** -2) - zeta ** 2) / zeta) - np.tan(
                2 * np.pi * f * zeta)

        # Initial guesses for zeta (we will use a range of guesses for multiple modes)
        zeta_initial_guesses = np.linspace(0.1, np.sqrt(H**2 * (beta1**-2 - beta2**-2)), 10)

        zeta_values = []
        if mode == 'newton':
            for zeta_guess in zeta_initial_guesses:
                try:
                    # Use root finding (newton-raphson) to find each mode (zeta)
                    zeta, _, _ = root_newton_raphson(zeta_guess, g, lambda zeta: (g(zeta + 1e-6) - g(zeta)) / 1e-6)
                    if zeta not in zeta_values:
                        zeta_values.append(zeta)  # Avoid duplicate zeta values
                except Exception:
                    continue  # Skip if root finding fails for a given initial guess

        elif mode == 'secant':
            for zeta_guess in zeta_initial_guesses:
                try:
                    # Use root finding (secant) to find each mode (zeta)
                    zeta, _, _ = root_secant_modified(zeta_guess, g, lambda zeta: (g(zeta + 1e-6) - g(zeta)) / 1e-6)
                    if zeta not in zeta_values:
                        zeta_values.append(zeta)  # Avoid duplicate zeta values
                except Exception:
                    continue  # Skip if root finding fails for a given initial guess

        # Compute c_L for each mode (zeta) and store it
        c_L_modes = [H / zeta for zeta in zeta_values]
        c_L_all_modes.append(c_L_modes)
        zeta_all_modes.append(zeta_values)

    return c_L_all_modes, zeta_all_modes

def mode_plt(c_L_all_modes, zeta_all_modes, frequencies):
    """ Plot Love wave velocity and wavelength for each mode. """
    plt.figure(figsize=(12, 8))

    # Plot c_L vs f for each mode
    plt.subplot(2, 1, 1)
    for i, c_L_modes in enumerate(c_L_all_modes):
        # Check if there are modes for this frequency
        if len(c_L_modes) > 0:
            plt.plot([frequencies[i]] * len(c_L_modes), c_L_modes, label=f"Frequency {frequencies[i]} Hz", marker='o', linestyle='None', color='blue')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Love wave velocity (m/s)')
    plt.title('Love Wave Velocity vs Frequency')
    plt.grid(True)
    plt.legend()

    # Plot wavelength (lambda_L) vs f for each mode
    plt.subplot(2, 1, 2)
    for i, c_L_modes in enumerate(c_L_all_modes):
        if len(c_L_modes) > 0:
            # Wavelength calculation: lambda_L = c_L / f
            lambda_L_modes = np.array(c_L_modes) / frequencies[i]
            plt.plot([frequencies[i]] * len(c_L_modes), lambda_L_modes, label=f"Frequency {frequencies[i]} Hz", marker='o', linestyle='None', color='red')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Wavelength (m)')
    plt.title('Wavelength vs Frequency')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


# Given parameters
rho1 = 1800  # kg/m^3
rho2 = 2500  # kg/m^3
beta1 = 1900  # m/s
beta2 = 3200  # m/s
H = 4000  # m (4 km)

# Optional: Plot asymptotes for root finding for some frequencies
atote_finding(frequencies, rho1, rho2, beta1, beta2, H)

# Find the dispersion curve and plot the results
c_L_all_modes, zeta_all_modes = dispersion_curve(frequencies, rho1, rho2, beta1, beta2, H)

# Plot the mode curves
mode_plt(c_L_all_modes, zeta_all_modes, frequencies)
