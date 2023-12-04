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
					["gain", ["gamma_mv", [1.1,0.2**2]], "continuous"], #mean, variance
					["rn",   ["gamma_mv", [2.5,0.25**2]], "continuous"], 
					["dc",   ["gamma_mv", [0.001,.001**2]], "continuous"]
				]
#need to update with range somehow? These can't be negative

fp_y_defs = ["y_gain", "y_rn", "y_dc"]

fp_d_defs = [
				["t_gain", ['uniform', [.1, 600]], "continuous"], #gain
				["I_gain", ['uniform', [1, 100]], "discrete"],    #gain
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
fp_x_defs = [
				#general
				["nx", ["nonrandom", [2048]], "discrete", 2048],
				["ny", ["nonrandom", [2048]], "discrete", 2048],
				["sigma_dc", ["uniform", [.3,.7]], "continuous", .5], #e-/s #WAGish based on a consistency test performed on SN20006
				["mu_stray", ["nonrandom", [0]], "continuous", 0], #e-/s #WAG for now
				["sigma_stray", ["uniform", [.001,.01]], "continuous", .005], #WAG for now
				#gain
				["P_signal", ["uniform", [.8,.95]], "continuous", 0.90], #Prob. of correctly identifying signal as event #WAG for now
				["P_noise", ["uniform", [.01,.1]], "continuous", 0.01], #Prob. of incorrectly identifying noise/interference as event #WAG for now
				["T_ccd", ["uniform", [_temp-1,_temp+1]], "continuous",  _temp], #K
				["E0", ["nonrandom", [_e0]], "continuous", _e0], #22.1 keV Cd-109 emission line
				["sigma_E", ["uniform", [math.sqrt((_m0 * _c**2) / (_k*(_temp+1)*_e0**2)),math.sqrt((_m0 * _c**2) / (_k*(_temp-1)*_e0**2))]], "continuous", _sigmaE], #1/eV^2
				["w", ["uniform", [3.64 + (2.12*.00025)*(300-_temp-1),3.70 + (2.80*.00025)*(300-_temp+1)]], "continuous", _w], #eV/e- #uncertainty in T, w(300K), a
				["activity_cd109", ["uniform", [1e-6,10e-6]], "continuous", 5e-6], #Ci #radioactivity of sample
				["grade_size", ["nonrandom", [3]], "discrete", 3], #3x3 event grade sizes
				["t_gain_setup", ["nonrandom", [1200]], "continuous", 1200], #WAG
				["t_gain_buffer", ["nonrandom", [5]], "continuous", 5], #WAG
				#rn
				["t_rn", ["nonrandom", [.1]], "continuous", .1], #100 ms exposure
				["t_rn_buffer", ["nonrandom", [5]], "continuous", 5], #WAG
				#dc
				["t_0", ["nonrandom", [.1]], "continuous", 0.1], #100ms baseline exposure assumed
				["t_dc_buffer", ["nonrandom", [5]], "continuous", 5], #WAG
				#qoi
				["tau", ["nonrandom", [1800]], "continuous", 1800],
				#cost
				["testbed_setup", ["nonrandom", [1800]], "continuous", 1800], #WAG
				#["C_engineer", ["nonrandom", [0.00694444444]], "continuous", 0.00694444444] #WAG $/s, from $25/hr
				["C_engineer", ["nonrandom", [1]], "continuous", 1] #just count time
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

	fp_d_defs.append([
						["n_qe", ['uniform', [0, 100]], "discrete"],   #qe
						["t_qe", ['uniform', [.1, 300]], "continuous"],#qe
						["I_qe", ['uniform', [1, 10]], "continuous"]   #qe  #WAG, check value
					])
		
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
					["S_pd", S_pd],   #representative photodiode, Hamamatsu S1337-1010BQ
					["S_pd_err", .01]  #mA/W
					#qoi
					#cost
					])


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
fp = ProblemDefinition(eta, H, Gamma, theta_req_defs, fp_y_defs, fp_d_defs, fp_x_defs)