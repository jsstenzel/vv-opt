import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv

sys.path.append('../..')
from problems.functionals import *
#llamas
from problems.snr_exp_models import *
from problems.snr_system_model import *

_wave_min = 350.0
_wave_max = 975.0
_bandpass = _wave_max - _wave_min

prior_gain = ["gamma_mv",  [1.1,0.2**2]] #mean, variance
prior_rn = ["gamma_mv", [2.5,0.25**2]]
prior_dc = ["gamma_mv", [0.001,.001**2]]
#the red ones might need a different prior than blue&green, based on the 2 test cameras
prior_gp_vph_red = ["gp_expquad", [prior_mean_vph_red, .0267, (_bandpass)**2]]  #mean_fn, variance, ls
prior_gp_vph_gre = ["gp_expquad", [prior_mean_vph_gre, .0267, (_bandpass)**2]]
prior_gp_vph_blu = ["gp_expquad", [prior_mean_vph_blu, .0267, (_bandpass)**2]]
prior_gp_sl = ["gp_expquad", [prior_mean_sl, .1, 1]]
prior_gp_bg = ["gp_expquad", [prior_mean_bg, .1, 1]]
prior_frd = ["gamma_mv", [0.077,0.022**2]]


#these priors are based on requirements that were met, see Camera Qual Report
theta_req_defs = [                             #mean, variance
					["gain_red", prior_gain, "continuous"],
					["gain_gre", prior_gain, "continuous"],
					["gain_blu", prior_gain, "continuous"],
					["rn_red", prior_rn, "continuous"],
					["rn_gre", prior_rn, "continuous"],
					["rn_blu", prior_rn, "continuous"],
					["dc_red", prior_dc, "continuous"],
					["dc_gre", prior_dc, "continuous"],
					["dc_blu", prior_dc, "continuous"],
					
					["vph_thru_red", prior_gp_vph_red, "continuous"],
					["vph_thru_gre", prior_gp_vph_gre, "continuous"],
					["vp_thruh_blu", prior_gp_vph_blu, "continuous"],
					["sl_thru_dichroic", prior_gp_sl, "continuous"],
					["bg_thru_dichroic", prior_gp_bg, "continuous"],
					["fiber_frd", prior_frd, "continuous"]
				]
#need to update with range somehow? These can't be negative

fp_y_defs = [	
				"y_gain_red", 
				"y_gain_gre", 
				"y_gain_blu", 
				"y_rn_red", 
				"y_rn_gre", 
				"y_rn_blu", 
				"y_dc_red", 
				"y_dc_gre", 
				"y_dc_blu", 
				
				"y_vph_red_pts", #expands?
				"y_vph_gre_pts", #expands?
				"y_vph_blu_pts", #expands?
				"y_sl_pts", #expands?
				"y_bg_pts", #expands?
				"y_frd"
			]

fp_d_defs = [
				["t_gain", ['uniform', [.1, 600]], "continuous"], #gain
				["I_gain", ['uniform', [1, 100]], "discrete"],    #gain
				["n_meas_rn", ['uniform', [1, 50]], "discrete"],  #rn
				["d_num", ['uniform', [2, 25]], "discrete"],      #dc
				["d_max", ['uniform', [1, 12000]], "continuous"], #dc
				["d_pow", ['uniform', [0,3]], "continuous"],      #dc
				
				["d_vph_n_pts", ['uniform', [0,_bandpass*2]], "discrete"],
				["d_dischroic_n_pts", ['uniform', [0,_bandpass*2]], "discrete"],
				["d_frd_n_meas", ['uniform', [0,2400]], "discrete"],
			]
	
_temp= -90+273.15 #K
_k = 1.380649e-23 #J / K
_c = 299792458 #m / s
_e0 = 22100 #eV
_m0 = 108.9049856 * 1.66054e-27 #cd109 atomic mass in kg
_sigmaE = math.sqrt((_m0 * _c**2) / (_k*_temp*_e0**2))
_w = 3.66 + 0.000615*(300-_temp)
fp_x_defs = [
				#focal plane
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
				#spectrograph
				
				#llamas
				
				#cost
				["testbed_setup", ["nonrandom", [1800]], "continuous", 1800], #WAG
				#["C_engineer", ["nonrandom", [0.00694444444]], "continuous", 0.00694444444] #WAG $/s, from $25/hr
				["C_engineer", ["nonrandom", [1]], "continuous", 1] #just count time
			]


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
llamas_snr = ProblemDefinition(eta, H, Gamma, theta_req_defs, fp_y_defs, fp_d_defs, fp_x_defs)


def update_llamas_problem(llamas_snr, d):
	d_masked = [(math.floor(dd) if self.d_masks[i]=='discrete' else dd) for i,dd in enumerate(d)]
	d_dict = dict(zip(self.d_names, d_masked))
	num_vph = d_dict["d_vph_n_pts"]
	num_dichroic = d_dict["d_dischroic_n_pts"]
	
	y_vph = ["y_vph_red_pts", "y_vph_gre_pts", "y_vph_blu_pts"]
	y_dichroic = ["y_sl_pts", "y_bg_pts"]
	
	mod_y_defs = []
	for yname in llamas_snr.y_names:
		if yname in y_vph:
			new_y = [yname+"_"+str(i) for i in range(num_vph)]
			mod_y_defs += new_y
		elif yname in y_dichroic: 
			new_y = [yname+"_"+str(i) for i in range(num_vph)]
			mod_y_defs += new_y
		else:
			mod_y_defs.append(yname)
	
	llamas_snr_d = ProblemDefinition(
		llamas_snr.eta, 
		llamas_snr.H, 
		llamas_snr.Gamma, 
		llamas_snr.theta_req_defs, 
		mod_y_defs, 
		llamas_snr.fp_d_defs, 
		llamas_snr.fp_x_defs
	)
		
	return llamas_snr_d