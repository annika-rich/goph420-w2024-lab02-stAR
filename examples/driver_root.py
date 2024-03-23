# Driver script for Lab02 implementation
from goph420_lab02.root_finding import root_newton_raphson
import numpy as np
import matplot.pyplot as plt

def main():
    """Main driver script for lab 2.
    """
    # constants of equation 1
    rho_1 = 1800 # kg / m ** 3
    rho_2 = 2500 # kg / m ** 3
    B1 = 1900 # m / s
    B2 = 3200 # m / s
    H = 4.0 # km

    # root finding equation g(zeta) = 0
    f = lambda f, zeta: (np.tan(2 * np.pi * f * zeta)
                        - rho_2 / rho_1 
                        * np.sqrt((H ** 2 * (B1 ** -2 - B2 ** -2) - zeta ** 2)) / zeta
    
                       )
    # derivative of root finding equation f
    dfdx = lambda f, zeta: 2 * np.pi * zeta * (1 / np.cos(2 * np.pi * zeta) ** 2)

if __name__ == "__main__":
    main()