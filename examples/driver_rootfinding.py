import numpy as np
import matplotlib.pyplot as plt
from src.goph420_lab03.roots import (root_secant_modified,
root_newton_raphson,
                                     )

# Frequency range
frequencies = [0.1, 0.5, 1, 1.4, 2, 5]  # Frequencies in Hz


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



# Dispersion equation function
def g(zeta, f, rho1, rho2, beta1, beta2, H):
    term1 = (rho2/rho1) * np.sqrt(H**2 * (1/beta1**2 - 1/beta2**2)-zeta**2)/zeta
    term2 = np.tan(2*np.pi*f*zeta)
    return term1 - term2

def func_deriv(zeta, f, rho1, rho2, beta1, beta2, H):
    term1 = (rho2/rho1)*(-1/np.sqrt((H**2*(beta1**-2 - beta2**-2))-zeta**2))
    term2 = ((rho2/rho1) * (-(np.sqrt(H**2 * (beta1**-2 - beta2**-2)-zeta**2))))/zeta**2
    term3 = 2*np.pi*f * (1/np.cos(2*np.pi*f*zeta))**2
    return term1 + term2 - term3

# Function to find cL (Love wave velocity) from zeta using Equation (2)
def compute_cL(zeta, beta1, H):
    cL_sq = 1/(beta1**-2 - (H**-2 / zeta**-2))
    return np.sqrt(cL_sq)


# Define a wrapper to pass additional parameters (f, rho1, rho2, beta1, beta2, H) to g and func_deriv
def wrapped_g(zeta, f, rho1, rho2, beta1, beta2, H):
    return g(zeta, f, rho1, rho2, beta1, beta2, H)

def wrapped_func_deriv(zeta, f, rho1, rho2, beta1, beta2, H):
    return func_deriv(zeta, f, rho1, rho2, beta1, beta2, H)


def find_all_modes(f, rho1, rho2, beta1, beta2, H, max_modes=5):

    zeta_max = np.sqrt(H ** 2 * (1 / beta1 ** 2 - 1 / beta2 ** 2))

    atotes = [0.0]
    k = 0
    while True:
        a = 0.25 * (2 * k + 1) / f
        if a < zeta_max:
            atotes.append(a)
            k += 1
        else:
            break
    atotes.append(zeta_max)

    zeta_roots = []

    for i in range(len(atotes) - 1):
        a, b = atotes[i], atotes[i + 1]
        zeta_guess = (b - 1e-4)
        print(zeta_guess)

        zeta_root, iter, _ = root_newton_raphson(
            zeta_guess,
            lambda zeta: g(zeta, f, rho1, rho2, beta1, beta2, H),
            lambda zeta: func_deriv(zeta, f, rho1, rho2, beta1, beta2, H)
        )

        if zeta_root > a and zeta_root < b and iter < 40:
            zeta_roots.append(zeta_root)


        if len(zeta_roots) >= max_modes:
            break

    return zeta_roots

# Modify the plot_love_wave_dispersion function to pass the required arguments
def plot_love_wave_dispersion_all_modes(frequencies, max_modes=5):
    mode_curves = [[] for _ in range(max_modes)]
    lambda_curves = [[] for _ in range(max_modes)]
    zeta_roots = []
    cL_roots = []
    lam_roots = []

    for f in frequencies:
        print(f"\n🔍 Solving for frequency: {f} Hz")

        zeta_new = find_all_modes(f, rho1, rho2, beta1, beta2, H, max_modes)
        cL_new = np.array([compute_cL(zeta, beta1, H) for zeta in zeta_new])
        lam_new = cL_new / f

        zeta_roots.append(zeta_new)
        cL_roots.append(cL_new)
        lam_roots.append(lam_new)

    for m in range(max_modes):
        for k,f in enumerate(frequencies):

            if m < len(cL_roots[k]):
                mode_curves[m].append(cL_roots[k][m])
                lambda_curves[m].append(lam_roots[k][m])



    print(mode_curves)
    print(lambda_curves)

    # Plot
    plt.figure(figsize=(14, 6))

# plot wave velocity
    plt.subplot(1, 2, 1)
    for mode_idx, cLs in enumerate(mode_curves):
        n_freq = len(cLs)
        plt.plot(frequencies[-n_freq: ], cLs, marker='o', label=f'Mode {mode_idx}')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Love Wave Velocity (m/s)")
    plt.title("Love Wave Velocity vs Frequency (All Modes)")
    plt.savefig("C:/Users/zcmit/git/Goph420projects/goph420-w2025-lab03-stZM/figures/wavevelocity-plot.png", dpi=300)
    plt.grid(True)
    plt.legend()

# plot wavelength
    plt.subplot(1, 2, 2)
    for mode_idx, lambdas in enumerate(lambda_curves):
        n_freq = len(lambdas)
        plt.plot(frequencies[-n_freq: ], lambdas, marker='x', label=f'Mode {mode_idx}')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Wavelength (m)")
    plt.title("Wavelength vs Frequency (All Modes)")
    plt.grid(True)
    plt.legend()

    plt.savefig("C:/Users/zcmit/git/Goph420projects/goph420-w2025-lab03-stZM/figures/wavelength-plot.png", dpi=300)
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

#plot cl vs f and lambda vs f
plot_love_wave_dispersion_all_modes(frequencies, max_modes=4)
