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

theta_req_defs = [ 
					["theta", ["uniform", [0,1]], "continuous"]
				]

y_defs = ["y"]

d_defs = [
				["d", ['uniform', [0, 1]], "continuous"], #gain
		 ]

x_defs = [
				#["nx", ["nonrandom", [2048]], "discrete", 2048],
			]

###Then, define the possible benchmarking models

def benchmark_1D_obs_model(theta, design, x, priormean=None, err=True):
	t = theta["theta"]
	d = design["d"]
	epsilon = 0
	if err:
		epsilon = scipy.stats.norm.rvs(scale = 10**(-4))
	y = (t**3)*(d**2) + t * math.exp(-abs(0.2-d)) + epsilon
	return [y]

def baseline_prediction_model(theta, x, verbose=False):
	#degenerate non-goal-oriented case
	return theta["theta"]

def zhongT1_prediction_model(theta, x, verbose=False):
	t = theta["theta"]
	z = math.sin(t) + t*math.exp(t + abs(0.5-t))
	return z
	
def zhongT2_prediction_model(theta, x, verbose=False):
	t = theta["theta"]
	if 0 <= t and t <= 0.15:
		return -100*t + 25
	elif 0.15 <= t and t <= 0.7:
		return 5
	elif 0.7 < t and t <= 1.0:
		return 50*t + 25
	else:
		return 0 #out of domain
		
def zhongT3_prediction_model(theta, x, verbose=False):
	t = theta["theta"]
	mu = 0.3
	sigma = 0.2
	z = (1/(math.sqrt(2*math.pi)*sigma)) * math.exp(-((t-mu)**2)/(2*(sigma**2)))
	return z
	
def bellcurve_prediction_model(theta, x, verbose=False):
	#for this problem, theta is the standard normal
	t = theta["theta"]
	#to get the desired distribution on Q, rescale theta and add nuisance uncertainty
	epsilon = scipy.stats.norm.rvs(scale = 5)
	
	q = 100 + t*10 + epsilon
	
	#None of this is actually what I want:
	#mu = .5
	#sigma = .025
	#inverse_cdf_theta = scipy.stats.norm.ppf(t, loc=mu, scale=sigma)
	#cdf_theta = scipy.stats.norm.cdf(t, loc=mu, scale=sigma)
	#pdf_theta = scipy.stats.norm.pdf(t, loc=mu, scale=sigma)
	
	return q
	
def longlefttail_prediction_mode(theta, x, verbose=False):
	#for this problem, theta is the standard normal
	t = theta["theta"]
	#now, nuisance uncertainty has long left tail
	epsilon = scipy.stats.norm.rvs(scale = 5) #TODO
	
	q = 100 + t*10 + epsilon
	return q

def longrighttail_prediction_model(theta, x, verbose=False):
	#for this problem, theta is the standard normal
	t = theta["theta"]
	#now, nuisance uncertainty has long right tail
	epsilon = scipy.stats.norm.rvs(scale = 5) #TODO
	
	q = 100 + t*10 + epsilon
	return q
	
def cost_simple(d, x):
	#penalize larger d linearly
	return 10*d[0]

#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
def make_1D_benchmark_problem(case=None):
	eta = benchmark_1D_obs_model
	Gamma = cost_simple
	tdefs = theta_req_defs
	ydefs = y_defs
	ddefs = d_defs
	xdefs = x_defs
	
	if case == "T0":
		H = baseline_prediction_model
	if case == "T1":
		H = zhongT1_prediction_model
	elif case == "T2":
		H = zhongT2_prediction_model
	elif case == "T3":
		H = zhongT3_prediction_model
	elif case == "H1":
		tdefs = [["theta", ["gaussian", [0,1]], "continuous"]]
		H = bellcurve_prediction_model
	else:
		H = baseline_prediction_model
		
	problem = ProblemDefinition(eta, H, Gamma, tdefs, ydefs, d_defs, x_defs)

	return problem