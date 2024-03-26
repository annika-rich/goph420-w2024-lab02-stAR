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
	tol = 1e-8 # approximate relative error tolerance
	maxit = 250 # max iterations
	eps_a = 2 * tol
	error = np.array([])
	itr = 1

	while eps_a > tol and itr < maxit:
		x1  = x0 - f(x0) / dfdx(x0) # Newton Raphson formula
		eps_a = np.abs((x1 - x0) / x1) # calculate approx. relative error
		error = np.append(error, eps_a)
		x0 = x1 # update guess
		itr += 1
		if itr >= maxit:
			raise RuntimeWarning(f"The Newton-Raphson Root Finding Algorithm did not converge over {itr} iterations.\nThe current value of x0: {x0}")
		
	return x0, itr, error

