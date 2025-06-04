#This details an obed script that can be called regardless of the problem
#TODO someday rework and generalize these methods to be independentof OBED approach, and even independent of the thing to be optimized

import sys
import numpy as np
import scipy.stats
import scipy.optimize
import matplotlib.pyplot as plt
from functools import partial
from tabulate import tabulate

#parallel:
from multiprocessing import Pool

sys.path.append('..')
from problems.problem_definition import *

#remove this later, not appropriate, this should be outside and passed in from a wrapper fn
from obed.obed_gbi import *

def parallel_design_print(costs, utils, designs, prob):
	headers = ["Cost", "Utility"]+list(prob.d_names)
	#print(res.history) #idk, returns None
	table = []
	ziplist = zip(costs, utils, designs)
	sorted_ziplist = sorted(ziplist, key=lambda x: x[1])
	for cost, util, design in sorted_ziplist:
		design_exact = prob._mask(design, prob.d_masks, prob.d_dists)
		row = [cost] + [util] + design_exact
		table.append(row)
	print(tabulate(table,headers),flush=True)

def minimize_with_penalty(problem, costcap, gmm_file, ylist_file, n_mc, n_tries, x0, ftol, penalty):
	print("Loading GMM and presamples...",flush=True)
	gmm = bn_load_gmm(gmm_file)
	presampled_ylist = bn_load_y(problem, ylist_file, doPrint=False, doDiagnostic=False)
	#for the initial guess, pick a pareto optimal design that satisfies the constraint
	d_lowers = [d_dist[1][0] for d_dist in problem.d_dists]
	d_uppers = [d_dist[1][1] for d_dist in problem.d_dists]
	bounds = scipy.optimize.Bounds(d_lowers, d_uppers)
	
	def fn_to_minimize(d):
		U_d,_ = U_varH_gbi_joint_presampled(d, problem, gmm, presampled_ylist, n_mc=n_mc, doPrint=False)
		C_d = problem.G(d)
		
		#Note that minimizing U_d actually maximizes the utility
		#Add a penalty factor if we exceed the cost cap
		if C_d <= costcap:
			return U_d
		else:
			return U_d + penalty + (costcap - C_d)
	
	costs = []
	utilities = []
	designs = []
	###do a few iterations of minimization here
	for i in range(n_tries):
		res = scipy.optimize.minimize(fn_to_minimize, x0, method='nelder-mead', options={'fatol': ftol, 'adaptive':True, 'disp': True}, bounds=bounds)
		costs.append(problem.G(res.x))
		utilities.append(res.fun)
		designs.append(res.x)
		
	parallel_design_print(costs, utilities, designs, problem)