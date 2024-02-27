import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv

sys.path.append('../..')
from problems.problem_definition import *
from problems.functionals import *

_wave_min = 350
_wave_max = 975
_bandpass = _wave_max - _wave_min

def prior_mean_vph_red(t):
	return 0

prior_gp_vph_red = ["gp_expquad", [prior_mean_vph_red, .05, (_bandpass)**2]]  #mean_fn, variance, ls


theta_defs = [
				["vph_thru_red", prior_gp_vph_red, "continuous"],
			]
#need to update with range somehow? These can't be negative

y_defs = 	[	
				"y_vph_red_pts", #expands?
			]

d_defs = 	[
				["d_vph_n_pts", ['uniform', [0,_bandpass*2]], "discrete"],
			]
	

x_defs = 	[
				["wave_min", ["nonrandom", [_wave_min]], "continuous", _wave_min],
				["wave_max", ["nonrandom", [_wave_max]], "continuous", _wave_max],
				["measurement_stddev", ["nonrandom", [.01]], "continuous", _wave_max]
			]


def eta(theta, d, x):
	vph_thru_red = theta["vph_thru_red"]
	d_vph_n_pts = d["d_vph_n_pts"]
	wave_min = x["wave_min"]
	wave_mix = x["wave_mix"]
	measurement_stddev = x["measurement_stddev"]
	
	#choose the measurement points
	measurement_pts = np.linspace(wave_min, wave_max, num=d_vph_n_pts)
	
	#make the measurements, assuming that there is one y_i for each measurement point ki
	y = []
	for ki in measurement_pts:
		yi = 0#do something with vph_thru_red
		y.append(yi)
		
	return y

def H(theta, x):
	vph_thru_red = theta["vph_thru_red"]
	
	#Take the average of theta over k:
	avg = 0
	
	return 0

def Gamma(d, x):
	d_vph_n_pts = theta["d_vph_n_pts"]
	return 100 + d_vph_n_pts #the more measurements, the more expensive


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
gp_test = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs, x_defs)


def update_gp_problem(problem, d):
	d_masked = [(math.floor(dd) if problem.d_masks[i]=='discrete' else dd) for i,dd in enumerate(d)]
	d_dict = dict(zip(problem.d_names, d_masked))
	num_vph = d_dict["d_vph_n_pts"]
	
	y_vph = ["y_vph_red_pts"]
	
	mod_y_defs = []
	for yname in problem.y_names:
		if yname in y_vph:
			new_y = [yname+"_"+str(i) for i in range(num_vph)]
			mod_y_defs += new_y
		else:
			mod_y_defs.append(yname)
	
	#problem_d = ProblemDefinition(
	#	problem.eta, 
	#	problem.H, 
	#	problem.G, 
	#	problem.theta_defs, 
	#	mod_y_defs, 
	#	problem.d_defs, 
	#	problem.x_defs
	#)
	
	problem_d = problem
	problem_d.dim_y = len(mod_y_defs)
	problem_d.y_names = mod_y_defs
		
	return problem_d
	
	
if __name__ == "__main__":
	print(gp_test)
	
	d = gp_test.sample_d(1)
	print('\n',d,'\n')
	
	llamas_snr_new = update_gp_problem(gp_test, d)
	
	print(llamas_snr_new)