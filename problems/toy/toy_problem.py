import sys
import os
import scipy.stats
import math

sys.path.append('../..')
#focal plane
from problems.problem_definition import *
from problems.fp_verification.fp_statistics import *
from problems.fp_verification.fp_experiment_models import *

useQE = False

def toy_eta(theta, d, x, err=True):
	t1 = theta["theta1"]
	d1 = d["d1"]
	random = scipy.stats.norm.rvs(loc=0, scale=d1, size=1)[0]
	return [t1 + random]
	
def toy_H(theta, x):
	H_scale = x["H_scale"]
	t1 = theta["theta1"]
	return H_scale*t1**2
	
def toy_Gamma(d, x):
	d1 = d["d1"]
	G_scale = x["G_scale"]
	return G_scale/d1
	
#def Ptheta_rvs(s=1):
#	#sample from the prior probability of theta
#	return scipy.stats.norm.rvs(size=s, loc=0.5, scale=0.5)
#
#def Ptheta_pdf(_theta):
#	return scipy.stats.norm.pdf(_theta, loc=0.5, scale=0.5)
#	
#def empty(whatever):
#	return 0

toy_theta_defs = [ 
					["theta1", ["gaussian", [0.5,0.5]], "continuous"], #mean, stddev
				 ]

toy_y_defs = ["y1"]

toy_d_defs = [
				["d1", ['uniform', [.05, .5]], "continuous"], #gain
			]
	
toy_x_defs = [
				#general
				["H_scale", ["uniform", [5,15]], "discrete", 10],
				["G_scale", ["uniform", [5,15]], "discrete", 10],
			]
			 


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
toy = ProblemDefinition(toy_eta, toy_H, toy_Gamma, toy_theta_defs, toy_y_defs, toy_d_defs, toy_x_defs)