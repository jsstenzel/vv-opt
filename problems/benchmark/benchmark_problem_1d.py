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

###Let's define a ProblemDefinition that tracks requirements
class ProblemReqDefinition(ProblemDefinition):
	def __init__(self, _eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs, _req=0, _req_type=None):
		super().__init__(_eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs)
		
		self.req = _req
		allowed_satisfaction_criteria = [
			"Q < r",
			"Q <= r",
			"Q > r",
			"Q >= r",
			"ra <= Q <= rb",
			"None"
		]
		if _req_type in allowed_satisfaction_criteria:
			self.req_type = _req_type
		elif _req_type == "ra <= Q <= rb" and not isinstance(_req, list):
			if len(_req) != 2:
				raise ValueError("For requirement satisfaction criterion ra <= Q <= rb, req must be a list of length 2. req: "+str(_req))
		else:
			raise ValueError("Unrecognized requirement satisfaction criterion: "+str(_req_type))

	def check_req(self, Q):
		if self.req_type == "Q < r":
			return Q < self.req
		elif self.req_type == "Q <= r":
			return Q <= self.req
		elif self.req_type == "Q > r":
			return Q > self.req
		elif self.req_type == "Q >= r":
			return Q >= self.req
		elif self.req_type == "ra <= Q <= rb":
			return (self.req[0] <= Q and Q <= self.req[1])
		else:
			return False
		

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
	
def zhongT1_gaussian(theta, x, verbose=False):
	z = zhongT1_prediction_model(theta, x, verbose=verbose)
	xi = scipy.stats.norm.rvs(scale = 0.1)
	return z + xi
	
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
	#T cases, from Zhong et al	
	if case == "T0":
		H = baseline_prediction_model
		eta = benchmark_1D_obs_model
		Gamma = cost_simple
		tdefs = [["theta", ["uniform", [0,1]], "continuous"]]
		ydefs = ["y"]
		ddefs = [["d", ['uniform', [0, 1]], "continuous"]]
		xdefs = []
		req = 0
		req_type = "None"
	elif case == "T1":
		H = zhongT1_prediction_model
		eta = benchmark_1D_obs_model
		Gamma = cost_simple
		tdefs = [["theta", ["uniform", [0,1]], "continuous"]]
		ydefs = ["y"]
		ddefs = [["d", ['uniform', [0, 1]], "continuous"]]
		xdefs = []
		req = 2.5
		req_type = "Q < r"
	elif case == "T2":
		H = zhongT2_prediction_model
		eta = benchmark_1D_obs_model
		Gamma = cost_simple
		tdefs = [["theta", ["uniform", [0,1]], "continuous"]]
		ydefs = ["y"]
		ddefs = [["d", ['uniform', [0, 1]], "continuous"]]
		xdefs = []
		req = 0
		req_type = "None"
	elif case == "T3":
		H = zhongT3_prediction_model
		eta = benchmark_1D_obs_model
		Gamma = cost_simple
		tdefs = [["theta", ["uniform", [0,1]], "continuous"]]
		ydefs = ["y"]
		ddefs = [["d", ['uniform', [0, 1]], "continuous"]]
		xdefs = []
		req = 0
		req_type = "None"

	#My own H cases - an adaptation of T1 with gaussian prior and noisy H
	elif case == "H1":
		req = 2.5
		req_type = "Q < r"
		H = zhongT1_gaussian
		eta = benchmark_1D_obs_model
		tdefs = [["theta", ["gaussian", [0.5,0.2]], "continuous"]]
		ddefs = [["d", ['uniform', [0, 1]], "continuous"]]
		ydefs = ["y"]
		xdefs = []
		Gamma = cost_simple
	elif case == "H2":
		#Same as H1, but with the requirement set to have 50/50 prior q
		req = 1.327
		req_type = "Q < r"
		H = zhongT1_gaussian
		eta = benchmark_1D_obs_model
		tdefs = [["theta", ["gaussian", [0.5,0.2]], "continuous"]]
		ddefs = [["d", ['uniform', [0, 1]], "continuous"]]
		ydefs = ["y"]
		xdefs = []
		Gamma = cost_simple
	elif case == "H3":
		#Same as H2, but with the requirement set to fail in the prior
		req = 0.75
		req_type = "Q < r"
		H = zhongT1_gaussian
		eta = benchmark_1D_obs_model
		tdefs = [["theta", ["gaussian", [0.5,0.2]], "continuous"]]
		ddefs = [["d", ['uniform', [0, 1]], "continuous"]]
		ydefs = ["y"]
		xdefs = []
		Gamma = cost_simple
	
	
	else:
		print("Don't recognize that case")
		sys.exit()
		
	problem = ProblemReqDefinition(eta, H, Gamma, tdefs, ydefs, ddefs, xdefs, req, req_type)

	return problem