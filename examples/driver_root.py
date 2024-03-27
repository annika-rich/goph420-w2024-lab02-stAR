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
    H = 4000 # m
    const = H ** 2 * ((B1 ** -2) - (B2 ** -2))
    # calculate the max possible value of zeta (for a real number result)
    zeta_max = np.sqrt(const)

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
        dx = -((2 * np.pi * f) / np.cos(2 * f * np.pi * zeta) ** 2
          + (rho_2 / rho_1) * np.sqrt(const - zeta ** 2) / zeta ** 2 
          + (rho_2 / rho_1) / (np.sqrt(const - zeta ** 2)))
        return dx

    modes = [0, 1, 2]
    # initiliaze empty 2D array to store root, velocity and wavelength at indices that indicate both the mode and frequency
    zeta = np.zeros([len(freq),len(modes)])
    c = np.zeros([len(freq), len(modes)])
    wvl = np.zeros([len(freq), len(modes)])
    m_max = []

    for i, f in enumerate(freq):
        plt.figure()
        plt.plot([0, zeta_max], [0, 0], '--b')
        z0 = 1e-4
        # initialize loop execution
        m = 0
        # initial guess calculated using locations of asymptotes of tan(2 * pi * f * zeta)
        # initial guess is below asymptote (below inflection point since using Newon_Raphson)
        func = lambda zeta: fx(f, zeta)
        deriv = lambda zeta: dx(f, zeta)
        x0 = (2 * m + 1) / (4 * f) - 1e-4
        while x0 < zeta_max:
            z1 = x0
            zp = np.linspace(z0, z1)
            #print(x0)
            plt.plot(zp, func(zp), '-k')
            plt.plot([z1 + 1e-4, z1 + 1e-4], [-100, 100], '--r')
            plt.plot(x0, 0.0, 'og')
            root, itr, error = root_newton_raphson(x0, func, deriv)
            plt.plot(root, func(root), 'ro')
            #print(root)
            if m < 3:
                zeta[i, m] = root
                c[i, m] = 1 / np.sqrt(B1 ** -2 - (zeta[i, m] / H ** 2))
                wvl[i, m] = c[i, m] / f
            m += 1
            # initialize new guess for next mode
            x0 = (2 * m + 1) / (4 * f) - 1e-4
            z0 = z1 + 2e-4
        
        if not m and func(x0:= zeta_max - 1e-4) < 0:
            root, itr, error = root_newton_raphson(x0, func, deriv)
            zeta[i, m] = root
            c[i, m] = 1 / np.sqrt(B1 ** -2 - (zeta[i, m] / H ** 2))
            wvl[i, m] = c[i, m] / f
            plt.plot(root, func(root), 'ro')
            m_max.append(m+1)
        else:
            m_max.append(m)

        z1 = zeta_max - 1e-4
        zp = np.linspace(z0, z1)
        plt.plot(zp, func(zp), '-k')
        plt.xlabel('zeta')
        plt.ylabel('f(zeta)')
        plt.ylim((-10, 10))
        plt.title(f"f = {f} Hz")
        plt.savefig(f'figures/frequency_{f}.png')

    print(c)
    print(wvl)
    print(np.argwhere(c[:,1] > 0).flatten())
        
    # make plots...

if __name__ == "__main__":
    main()