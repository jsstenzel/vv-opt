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

useQE = False

eta = fp_likelihood_fn
H = fp_hlva
Gamma = fp_cost_simple			   

#these priors are based on requirements that were met, see Camera Qual Report
theta_req_defs = [ 
					("gain", ["uniform", [1.0-0.2,1.0+0.2]]), #mean, stddev
					("rn",   ["uniform", [2.5-0.25,2.5+0.25]]), #these should probably be gamma fns
					("dc",   ["uniform", [0.001-.0002,0.001+.003]]),
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
				("mu_stray", .1), #e-/s #WAG for now
				("sigma_stray", .005), #WAG for now
				#gain
				("P_signal", 0.95), #Prob. of correctly identifying signal as event #WAG for now
				("P_noise", 0.01), #Prob. of incorrectly identifying noise/interference as event #WAG for now
				("T_ccd", _temp), #K
				("E0", _e0), #22.1 keV Cd-109 emission line
				("sigma_E", math.sqrt((_m0 * _c**2) / (_k*_temp*_e0**2))), #1/eV^2
				("w", 3.66 + 0.000615*(300-_temp)), #eV/e- #this has uncertainty. nuisance parameter?
				("rate_cd109", 1.7387e-8), #1/s #based on 461.4 days Cd-109 half life
				("grade_size", 3), #3x3 event grade sizes
				("t_gain_setup", 1200), #WAG
				("t_gain_buffer", 5), #WAG
				#rn
				("t_rn", .1), #100 ms exposure
				("t_rn_buffer", 5), #WAG
				#dc
				("t_0", 0.1), #100ms baseline exposure assumed
				("t_dc_buffer", 5), #WAG
				#qoi
				("tau", 1800),
				#cost
				("testbed_setup", 1800), #WAG
				("C_engineer", 0.00694444444) #WAG $/s, from $25/hr
			 ]
			 
if useQE == True:
	eta = fp_qe_likelihood_fn
	H = fp_qe_hlva
	Gamma = fp_qe_cost			   

	#this prior is based on measurements for red CCD42-40
	theta_req_defs.append(
		("qe", ["funct_splines", [
					[(400,.25),(500,.45),(650,.75),(900,.45),(975,.05)], #measured data
					3, #order of interpolation
					.05, #y error on all points
					[350,975], #LLAMAS spectral range
					[0.0,1.0], #range of QE
			 ]]
		))

	fp_y_defs.append("y_qe")

	fp_d_defs.extend(["n_qe", "t_qe", "I_qe", "S_err"]) #qe design variables
		
	S_pd = Functional([(200,.12),(280,.1),(300,.125),(400,.185),(633,.33),(930,.5),(1000,.45),(1100,.15)]) #nm, A/W
	S_pd.set_xlim(350,975)
	S_pd.set_ylim(0,1)
	S_pd.spline_interp(3)
	fp_x_defs.extend([
					#general
					#gain
					#rn
					#dc
					#qe
					("S_pd", S_pd),   #representative photodiode, Hamamatsu S1337-1010BQ
					("S_pd_err", .01)  #mA/W
					#qoi
					#cost
				    ])


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
fp = ProblemDefinition(eta, H, Gamma, theta_req_defs, fp_y_defs, fp_d_defs, fp_x_defs)

#print(fp)
#test1 = fp.prior_rvs(1)
#test3 = fp.prior_rvs(3)
#print(test1)
#print(fp.prior_pdf_unnorm(test1))
#print(test3)