import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
import csv

sys.path.append('../..')
#focal plane
from problems.problem_definition import *
from problems.fp_verification.fp_statistics import *
from problems.fp_verification.fp_experiment_models import *

def fp_hlva_2(theta, x, verbose=False):
	if type(theta) is dict:
		new_theta = [0, theta["rn"], theta["dc"]]
	else:
		new_theta = [0, theta[0], theta[1]]
	return fp_hlva(new_theta, x, verbose=verbose)

def fp_likelihood_2(theta, d, x, dummy, err=True):
	#define interest params:
	if type(theta) is dict:
		rn = theta["rn"]
		dc = theta["dc"]
	else:
		rn = theta[0]
		dc = theta[1]
	#define design variables, enforcing discreteness:
	n_meas_rn = d["n_meas_rn"]
	d_num = d["d_num"]
	d_max = d["d_max"]
	d_pow = d["d_pow"]
	#just pass along entire x
	gain = x["gain"]
	
	y2 = read_noise_exp(gain, rn, n_meas_rn, x, err)
	y3 = dark_current_exp(gain, rn, dc, d_num, d_max, d_pow, x, err)
	
	return [y2, y3]
 
def fp_cost_2(d, x):
	d_old = d
	d_old["t_gain"] = 235.632 #dstar
	d_old["I_gain"] = 92 #dstar
	return fp_cost_simple(d, x)

#these priors are based on requirements that were met, see Camera Qual Report
fp2_theta_defs = [ 
					["rn",   ["gamma_mv", [2.5,0.25**2]], "continuous"], 
					["dc",   ["gamma_mv", [0.001,.001**2]], "continuous"]
				]

fp2_y_defs = ["y_rn", "y_dc"]

fp2_d_defs = [
				["n_meas_rn", ['uniform', [1, 50]], "discrete"],  #rn
				["d_num", ['uniform', [2, 25]], "discrete"],     #dc
				["d_max", ['uniform', [1, 12000]], "continuous"], #dc
				["d_pow", ['uniform', [0,3]], "continuous"]       #dc
			]
	
_temp= -90+273.15 #K
_k = 1.380649e-23 #J / K
_c = 299792458 #m / s
_e0 = 22100 #eV
_m0 = 108.9049856 * 1.66054e-27 #cd109 atomic mass in kg
_sigmaE = math.sqrt((_m0 * _c**2) / (_k*_temp*_e0**2))
_w = 3.66 + 0.000615*(300-_temp)
fp2_x_defs = [
				#sequential
				["gain", ["gamma_mv", [1.0999913287843317,1.5998220769876884e-06]], "continuous", 1.1], #gain posterior mean and variance
				#general
				["nx", [], "discrete", 2048],
				["ny", [], "discrete", 2048],
				["sigma_dc", [], "continuous", .5], #e-/s #WAGish based on a consistency test performed on SN20006
				["mu_stray", [], "continuous", 0], #e-/s #WAG for now
				["sigma_stray", [], "continuous", .005], #WAG for now
				#gain
				["P_signal", [], "continuous", 0.90], #Prob. of correctly identifying signal as event #WAG for now
				["P_noise", [], "continuous", 0.01], #Prob. of incorrectly identifying noise/interference as event #WAG for now
				["T_ccd", [], "continuous",  _temp], #K
				["E0", [], "continuous", _e0], #22.1 keV Cd-109 emission line
				["sigma_E", [], "continuous", _sigmaE], #1/eV^2
				["w", [], "continuous", _w], #eV/e- #uncertainty in T, w(300K), a
				["activity_cd109", [], "continuous", 5e-6], #Ci #radioactivity of sample
				["grade_size", [], "discrete", 3], #3x3 event grade sizes
				["t_gain_setup", [], "continuous", 1200], #WAG
				["t_gain_buffer", [], "continuous", 5], #WAG
				#rn
				["t_rn", [], "continuous", .1], #100 ms exposure
				["t_rn_buffer", [], "continuous", 5], #WAG
				#dc
				["t_0", [], "continuous", 0.1], #100ms baseline exposure assumed
				["t_dc_buffer", [], "continuous", 5], #WAG
				#qoi
				["tau", [], "continuous", 1800],
				#cost
				["testbed_setup", [], "continuous", 1800], #WAG
				["C_engineer", [], "continuous", 1] #just count time
			]

#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
fp2 = ProblemDefinition(fp_likelihood_2, fp_hlva_2, fp_cost_2, fp2_theta_defs, fp2_y_defs, fp2_d_defs, fp2_x_defs)

if __name__ == '__main__':  
	print(fp2)
	print(fp2.sample_x(5))
