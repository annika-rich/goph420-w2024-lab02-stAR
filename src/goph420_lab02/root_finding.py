# Utility functions for Lab02
import numpy as np 

def root_newton_raphson(x0, f, dfdx):
	"""This function implements the Newton-Raphson Root Finding Open Method.

	Parameters
	----------
	x0:     float
	initial guess for value of the root
	f:      function
	Function for which the root is evaluated
	dfdx:   function
	first derivative of function, f, for which the root is evaluated

	Returns
	-------
	float
	final estimate of the root of the given function

	int
	number of iterations to convergence

	numpy.ndarray
	1D vector of approximate relative error at each iteration
	"""
	# initialize loop variables
	tol = 1e-16 # approximate relative error toleranc
	maxit = 100 # max iterations
	esp_a = 2 * tol
	itr = 1

	while eps_a > tol & itr < maxit:
		x1  = x0 - f(x0) / dfdx(x0) # Newton Raphson formula
		eps_a = np.abs((x1 - x0) / x1) # calculate approx. relative error
		x0 = x1 # update guess
		if itr >= maxit:
			raise RuntimeWarning(f"The Newton-Raphson Root Finding Algorithm did not converge over {itr} iterations.")
		


