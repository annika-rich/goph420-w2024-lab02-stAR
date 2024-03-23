import unittest
from scipy.optimize import newton
import numpy as np

from goph420_lab02.root_finding import root_newton_raphson


class TestNewtonRaphson(unittest.TestCase):

    def setUp(self):
        bracket = -0.075
        r = 1.75
        self.h = 0.6
        self.f = lambda h: 4 * bracket * r ** 3 + 3 * r * h ** 2 - h ** 3
        self.dfdx = lambda h : -3 * (h - 3.5) * h
    
    # from Quiz 2 
    def test_value(self):
        excel_root = 0.5872119182
        excel_itr = 5
        excel_error = np.array([2.15888936456e-2, 1.84730809558e-4, 1.36220772173e-8, 1.89066841166e-16])
        root, itr, error = root_newton_raphson(self.h, self.f, self.dfdx)
        self.assertAlmostEqual(root, excel_root)
        self.assertAlmostEqual(itr, excel_itr)
        self.assertTrue(np.allclose(error, excel_error, atol = 1e-8))
