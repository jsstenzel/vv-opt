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

fp_theta_defs = [ 
					("gain", ["gaussian", [0,1]]),
					("rn",   ["gaussian", [0,1]]),
					("dc",   ["gaussian", [0,1]]),
				]

fp_y_defs = ["y_gain", "y_rn", "y_dc"]

fp_d_defs = [
				"t_gain", "I_gain",        #gain
				"n_meas_rn",               #rn
				"d_num", "d_max", "d_pow"  #dc
			 ]
	
fp_x_defs = [
				#general
				("nx", 2048),
				("ny", 2048),
				("sigma_dc", .5), #e-/s
				("sigma_stray", .1) ,
				#gain
				("P_signal", *),
				("P_noise", *),
				("T_ccd", -90+273.15), #K
				("E0", *),
				("sigma_E", *),
				("w", *), #i probably need to treat this as a nuisance parameter later?
				("rate_cd109", *), #look at cd-109 sample. assume 1 gram maybe?
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
fp = ProblemDefinition(fp_likelihood_fn, fp_hlva, fp_cost, fp_theta_defs, fp_y_defs, fp_d_defs, fp_x_defs)

print(fp)