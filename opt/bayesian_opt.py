#This details an optimizer than can be called regardless of the problem

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial

from skopt.plots import plot_gaussian_process
from skopt import gp_minimize
from skopt.plots import plot_convergence

sys.path.append('..')
from problems.problem_definition import *

#here, f is 
#

"""
problem			The ProblemDefinition object
f				a function that takes a single vector and returns a single scalar
				for design optimization, the input may be a design and the output may be a utility
n_init			number of points to evaluate f(x) before running optimizer
				you want this high enough to explore the design space
n_iter			number of steps you take towards global optimum
f_noise			if you've done uncertainty analysis for f, use the variance of f here
findMinimum		set to False to maximize f(x) instead
plotConv		set to True to see convergence plot
"""
def bayesian_opt(problem, f, n_init, n_iter, f_noise="gaussian", findMinimum=True, plotConv=False):

	d_lowers = [d_dist[1][0] for d_dist in problem.d_dists]
	d_uppers = [d_dist[1][1] for d_dist in problem.d_dists]
	d_bounds = [(d_lower,d_upper) for d_upper,d_lower in zip(d_uppers,d_lowers)]
	
	def negative_f(x):
		return -f(x)
	
	if findMinimum: #handle both minimizing and maximizing
		opt_f = f
	else:
		opt_f = negative_f
						
	print("Finding",("minimum" if findMinimum else "maximum"),"of",f,"...",flush=True)
	np.int = int
	res = gp_minimize(opt_f,			  # the function to minimize
					d_bounds,		   # the bounds on each dimension of x
					acq_func="EI",	  # the acquisition function
					n_calls=n_iter,	# the number of evaluations of f
					n_initial_points=n_init,  # the number of random initialization points
					initial_point_generator="lhs", 
					noise=f_noise, 			#0.1**2,	   # the noise level -- 
					verbose=True)
					#random_state=1234)  # the random seed
	
	print(res)
	optimal_d = res.x
	best_f = res.fun if findMinimum else -res.fun #handle both minimizing and maximizing
	
	if plotConv:
		plot_convergence(res)
		
	return optimal_d, best_f
	