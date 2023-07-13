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
from problems.fp_verification.fp_statistics import *
from problems.fp_verification.fp_experiment_models import *


#eta = fp_likelihood_fn
#H = fp_hlva
#Gamma = fp_cost			   

#these priors are based on requirements that were met, see Camera Qual Report
theta_req_defs = [ 
					("gain", ["gaussian", [1.0,0.2]]), #mean, stddev
					("rn",   ["gaussian", [2.5,0.25]]),
					("dc",   ["gaussian", [0.001,.002]]),
				]
#need to update with range somehow? These can't be negative

#for post-fabrication pre-testing priors, I should update the above priors
#with the Qual data from SN20001 and SN20003
#theta_qual_defs = [ 
#					("gain", ["gaussian", [1.0,0.2]]), #mean, stddev
#					("rn",   ["gaussian", [2.5,0.25]]),
#					("dc",   ["gaussian", [0.001,.002]]),
#				]

fp_y_defs = ["y_gain", "y_rn", "y_dc"]

fp_d_defs = [
				"t_gain", "I_gain",        #gain
				"n_meas_rn",               #rn
				"d_num", "d_max", "d_pow"  #dc
			 ]
	
_temp= -90+273.15 #K
_k = 1.380649e-23 #J / K
_c = 299792458 #m / s
_e0 = 22100 #eV
_m0 = 108.9049856 * 1.66054e-27 #cd109 atomic mass in kg
fp_x_defs = [
				#general
				("nx", 2048),
				("ny", 2048),
				("sigma_dc", .5), #e-/s #WAG for now
				("sigma_stray", .1), #WAG for now
				#gain
				("P_signal", 0.95), #Prob. of correctly identifying signal as event #WAG for now
				("P_noise", 0.01), #Prob. of incorrectly identifying noise/interference as event #WAG for now
				("T_ccd", _temp), #K
				("E0", e0), #22.1 keV Cd-109 emission line
				("sigma_E", math.sqrt((_m0 * _c^2) / (_k*_temp*_e0^2))), #1/eV^2
				("w", 3.66 + 0.000615*(300-_temp)), #eV/e- #this has uncertainty. nuisance parameter?
				("rate_cd109", 1.7387e-8), #1/s #based on 461.4 days Cd-109 half life
				("grade_size", 3), #3x3 event grade sizes
				#rn
				("t_rn", .1), #100 ms exposure
				#dc
				("t_0", 0.1), #100ms baseline exposure assumed
				#qoi
				("tau", 1800),
				#cost
			 ]


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
fp = ProblemDefinition(fp_likelihood_fn, fp_hlva, fp_cost, theta_req_defs, fp_y_defs, fp_d_defs, fp_x_defs)

print(fp)