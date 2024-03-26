# Driver script for Lab02 implementation
from goph420_lab02.root_finding import root_newton_raphson
import numpy as np
import matplotlib.pyplot as plt

def main():
    """Main driver script for lab 2.
    """
    # constants of equation 1
    rho_1 = 1800 # kg / m ** 3
    rho_2 = 2500 # kg / m ** 3
    B1 = 1900 # m / s
    B2 = 3200 # m / s
    H = 4.0 # km
    const = H ** 2 * ((B1 ** -2) - (B2 ** -2))
    # calculate the max possible value of zeta (for a real number result)
    zeta_max = np.sqrt(const)
    print(f"zeta max: {zeta_max}")
    # define list of possible love wave frequencies to plot
    freq = np.array([0.01, 0.1, 0.5, 1, 5, 10]) # Hz
     # root finding equation g(zeta) = 0
    def fx(f, zeta):
        fx = (rho_2 / rho_1 
         * np.sqrt((const - zeta ** 2)) / zeta
         - np.tan(2 * np.pi * f * zeta))
        return fx
        
    # derivative of root finding equation f
    def dx(f, zeta):
        dx = ((2 * np.pi * f) / np.cos(2 * np.pi * zeta)
          + (rho_2 / rho_1) * np.sqrt(const - zeta ** 2) / zeta ** 2 
          + (rho_2 / rho_1) * 1 / (np.cos(const - zeta ** 2)))
        return dx

    modes = np.array([0, 1, 2])
    # initiliaze empty 2D array to store root at indices that indicate both the mode and frequency
    zeta = np.zeros([len(freq),len(modes)])

    for i, f in enumerate(freq):
        for m in modes:
            # initial guess calculated using locations of asymptotes of tan(2 * pi * f * zeta)
            # initial guess is below asymptote (below inflection point since using Newon_Raphson)
            x0 = (2 * m + 1) / (4 * f) - 1e-4
            func = lambda zeta: fx(f, zeta)
            deriv = lambda zeta: dx(f, zeta)
            root, itr, error = root_newton_raphson(x0, func, deriv)
            zeta[i, m] = root
            


if __name__ == "__main__":
    main()