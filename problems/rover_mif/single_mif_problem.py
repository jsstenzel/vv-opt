import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
import csv

sys.path.append('../..')
from problems.problem_definition import *

###First, define the fixed def objects

#probability of a mission-interrupting fault
theta_defs = [ 
					["motor drive electronics fail rate [per hr]", ["gamma_ab", [2,800]], "continuous"] #TODO values
				]

y_defs = [
					"n_faults_detected"
		]

d_defs = [
					["rover test time [hr]", ['uniform', [0, 10000]], "continuous"], #TODO values
		 ]

x_defs = []

###Then, define the possible benchmarking models

#note that y is discrete, result of a Poisson process
def rover_drive_test(theta, design, x, priormean=None, err=True):
	t = theta["motor drive electronics fail rate [per hr]"]
	d = design["rover test time [hr]"]
	
	#sample from Poisson process
	y = np.random.poisson(t*d, 1)
	
	return y

def H_degenerate(theta, x, verbose=False):
	t = theta["motor drive electronics fail rate [per hr]"]
	return t

def cost_simple(d):
	return d
	
class ProblemDefinitionPoissonGaussian(ProblemDefinition):
	def __init__(self, _eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs):
		super().__init__(_eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs)
		
	def sample_posterior(self, y, d):
		vals = [] #a list length num_vals of random numbers of size dim_theta
		for prior in self.priors: ###iterate over dim_theta
			dtype = prior[0]
			params = prior[1]
			#need to do this carefully, we have multiple thetas and multiple samples
	
			if dtype == 'gamma_ab':
				alpha = params[0] + y
				beta = params[1] + d
				post_thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha, scale=1.0/beta)
				vals.append(post_thetas_i.tolist())

		#turn from list of rvs at each prior, to list of theta rvs
		vals = np.transpose(vals)
	
		#return
		if len(vals) == 1:
			return vals[0] #this is a list (dim_theta)
		else:
			return vals #this is a list (num_vals) of random variables (dim_theta)

	def singlemif_infocriterion(self, d, y):
		params = self.priors[0][1]
		a = params[0]
		b = params[1]
		u = a*np.log((b+d)/b) - scipy.special.gammaln(a+y) + scipy.special.gammaln(a) + y*scipy.special.digamma(a+y) - d*((a+y)/(b+d))
		return u
		
	def utility(self, d, n_MC):
		thetas = self.prior_rvs(n_MC)
		ys = [self.eta([theta], [d]) for theta in thetas]
		U = 0
		for i,y in enumerate(ys):
			ui = self.singlemif_infocriterion(d, y[0]) #different utility for each theta,y
			U += ui
		return U/n_MC

#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
def make_single_mif_problem():
	eta = rover_drive_test
	Gamma = cost_simple
	tdefs = theta_defs
	ydefs = y_defs
	ddefs = d_defs
	xdefs = x_defs
	H = H_degenerate
		
	problem = ProblemDefinitionPoissonGaussian(eta, H, Gamma, tdefs, ydefs, d_defs, x_defs)

	return problem