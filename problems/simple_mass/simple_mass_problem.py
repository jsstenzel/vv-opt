import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv

sys.path.append('../..')
#focal plane
from problems.problem_definition import *


#Define the mass problem
#But allow some hyperparameters!
#alphas are the prior mean masses, betas are the prior stddev masses
def simple_mass_problem_def(alphas, betas):
	if len(alphas) != len(betas):
		print("simple_mass_problem_def error, alphas and betas dont match")
		sys.exit()
	N = len(alphas)
	
	#Each parameter of interest is the mass of an element in the system
	theta_dists = [ ["gaussian", [alphas[i],betas[i]]] for i in range(N)]
	theta_defs = [ ["t"+str(i),theta_dists[i],"continuous"] for i in range(N)]

	y_defs = ["y"+str(i) for i in range(N)]

	#Each design variable is a measurement precision
	d_dist = ["uniform", [0.001,10]] #i suppose?
	d_defs = [ ["d"+str(i),d_dist,"continuous"] for i in range(N)]

	x_defs = [
		["meas_err_scale", [], "continuous", 1.0],
		["cost_scale", [], "continuous", 1.0],
	]

	def eta(_t, _d, _x, err=True):
		if err:
			err_scale = _x["meas_err_scale"]
		else:
			err_scale = 0
		#generalize? ugh i wish this wasnt a dictionary
		ds = [ _d["d"+str(i)] for i in range(N)]
		ts = [ _t["t"+str(i)] for i in range(N)]
		
		y = []
		for i,t in enumerate(ts):
			epsilon = scipy.stats.norm.rvs(scale = ds[i] * err_scale)
			yi = t + ds[i] * epsilon
			y.append(yi)
		return y
		
	#Goal is to determine accurate measurement of total mass
	def H(_t, _x):
		ts = [ _t["t"+str(i)] for i in range(N)]
		total_mass = sum([ _t["t"+str(i)] for i in range(N)])
		
		return total_mass
		
	#Cost scales inversely with experimental precision
	def Gamma(_d, _x):
		ds = [ _d["d"+str(i)] for i in range(N)]
		costs = [_x["cost_scale"]/d for d in ds]
		total_cost = sum(costs)
		
		return total_cost		   


	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	simple_mass = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs, x_defs)
	return simple_mass