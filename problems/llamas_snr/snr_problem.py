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
from problems.gaussian_process import *
#llamas
from problems.llamas_snr.snr_exp_models import *
from problems.llamas_snr.snr_system_model import *

#set path variables
_dir = "./llamas-etc/COATINGS/"

_wave_min = 350.0
_wave_max = 975.0
_bandpass = _wave_max - _wave_min
#my physical argument here is that we should consider measurements at these distances apart to be essentially independent
#and this is based on the approximate spectral resolution of the instrument
#R = lambda / dlambda = 2200
#so for lambda=350, dlambda=0.159
#and for lambda=975, dlambda=0.443
_lengthscale = 1.0

prior_gain = ["gamma_mv",  [1.1,0.2**2]] #mean, variance
prior_rn = ["gamma_mv", [2.5,0.25**2]]
prior_dc = ["gamma_mv", [0.001,.001**2]]
#the red ones might need a different prior than blue&green, based on the 2 test cameras
ppts, meanfn = get_ppts_meanfn_file(_dir+"CCD42-40_dd.txt", 3)
prior_qe_red = ["gp_expquad", [.05, _lengthscale, ppts, meanfn]]
ppts, meanfn = get_ppts_meanfn_file(_dir+"CCD42-40_green.txt", 3)
prior_qe_gre = ["gp_expquad", [.05, _lengthscale, ppts, meanfn]]
ppts, meanfn = get_ppts_meanfn_file(_dir+"CCD42-40_blue.txt", 3)
prior_qe_blu = ["gp_expquad", [.05, _lengthscale, ppts, meanfn]]
	
ppts, meanfn = get_ppts_meanfn_file(_dir+"wasach_llamas2200_red.txt", 2)
prior_gp_vph_red = ["gp_expquad", [.0267, _lengthscale, ppts, meanfn]]
ppts, meanfn = get_ppts_meanfn_file(_dir+"wasach_llamas2200_green.txt", 2)
prior_gp_vph_gre = ["gp_expquad", [.0267, _lengthscale, ppts, meanfn]]
ppts, meanfn = get_ppts_meanfn_file(_dir+"wasach_llamas2200_blue.txt", 2)
prior_gp_vph_blu = ["gp_expquad", [.0267, _lengthscale, ppts, meanfn]]

ppts, meanfn = get_ppts_meanfn_file(_dir+"ECI_FusedSilica.txt", 3)
prior_gp_sl = ["gp_expquad", [.1, _lengthscale, ppts, meanfn]]
ppts, meanfn = get_ppts_meanfn_file(_dir+"ECI_FusedSilica.txt", 3)
prior_gp_bg = ["gp_expquad", [.1, _lengthscale, ppts, meanfn]]

prior_frd = ["gamma_mv", [0.077,0.022**2]]


#these priors are based on requirements that were met, see Camera Qual Report
theta_defs = [                             #mean, variance
					["gain_red", prior_gain, "continuous"],
					["gain_gre", prior_gain, "continuous"],
					["gain_blu", prior_gain, "continuous"],
					["rn_red", prior_rn, "continuous"],
					["rn_gre", prior_rn, "continuous"],
					["rn_blu", prior_rn, "continuous"],
					["dc_red", prior_dc, "continuous"],
					["dc_gre", prior_dc, "continuous"],
					["dc_blu", prior_dc, "continuous"],
					["qe_red", prior_qe_red, "functional"],
					["qe_gre", prior_qe_gre, "functional"],
					["qe_blu", prior_qe_blu, "functional"],
					
					["vph_thru_red", prior_gp_vph_red, "functional"],
					["vph_thru_gre", prior_gp_vph_gre, "functional"],
					["vph_thru_blu", prior_gp_vph_blu, "functional"],
					["sl_thru_dichroic", prior_gp_sl, "functional"],
					["bg_thru_dichroic", prior_gp_bg, "functional"],
					["fiber_frd", prior_frd, "continuous"]
				]
#need to update with range somehow? These can't be negative

y_defs = [	
				"y_gain_red", 
				"y_gain_gre", 
				"y_gain_blu", 
				"y_rn_red", 
				"y_rn_gre", 
				"y_rn_blu", 
				"y_dc_red", 
				"y_dc_gre", 
				"y_dc_blu", 
				"y_qe_red_pts", #expands
				"y_qe_gre_pts", #expands
				"y_qe_blu_pts", #expands
				
				"y_vph_red_pts", #expands
				"y_vph_gre_pts", #expands
				"y_vph_blu_pts", #expands
				"y_sl_pts", #expands
				"y_bg_pts", #expands
				"y_frd"
			]

d_defs = [
				["t_gain", ['uniform', [.1, 600]], "continuous"], #gain
				["I_gain", ['uniform', [1, 100]], "discrete"],    #gain
				["n_meas_rn", ['uniform', [1, 50]], "discrete"],  #rn
				["d_num", ['uniform', [2, 25]], "discrete"],      #dc
				["d_max", ['uniform', [1, 12000]], "continuous"], #dc
				["d_pow", ['uniform', [0,3]], "continuous"],      #dc
				
				["n_qe", ['uniform', [0, 100]], "discrete"],   #qe
				["t_qe", ['uniform', [.1, 600]], "continuous"],#qe
				["I_qe", ['uniform', [1, 10]], "continuous"],   #qe  #WAG, check value
				
				["d_vph_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
				["d_dichroic_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
				["d_frd_n_meas", ['uniform', [0,2400]], "discrete"],
			]
	
_temp= -90+273.15 #K
_k = 1.380649e-23 #J / K
_c = 299792458 #m / s
_e0 = 22100 #eV
_m0 = 108.9049856 * 1.66054e-27 #cd109 atomic mass in kg
_sigmaE = math.sqrt((_m0 * _c**2) / (_k*_temp*_e0**2))
_w = 3.66 + 0.000615*(300-_temp)
photodiode = [(200,.12),(280,.1),(300,.125),(400,.185),(633,.33),(930,.5),(1000,.45),(1100,.15)] #representative photodiode, Hamamatsu S1337-1010BQ

#build model
spectrograph_models = build_model(_dir)
skyfile_path = os.environ['COATINGS_PATH']+"eso_newmoon_radiance.txt"

x_defs = [
				#focal plane
				["nx", ["nonrandom", [2048]], "discrete", 2048],
				["ny", ["nonrandom", [2048]], "discrete", 2048],
				["pixel_size_mm", [], "continuous", 13.5],
				["sigma_dc", ["uniform", [.3,.7]], "continuous", .5], #e-/s #WAGish based on a consistency test performed on SN20006
				["mu_stray", ["nonrandom", [0]], "continuous", 0], #e-/s #WAG for now
				["sigma_stray", ["uniform", [.001,.01]], "continuous", .005], #WAG for now
				#gain exp
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
				#rn exp
				["t_rn", [], "continuous", .1], #100 ms exposure
				["t_rn_buffer", [], "continuous", 5], #WAG
				#dc exp
				["t_0", [], "continuous", 0.1], #100ms baseline exposure assumed
				["t_dc_buffer", [], "continuous", 5], #WAG
				#qe
				["S_pd", [], "functional", define_functional([p[0] for p in photodiode], [p[1] for p in photodiode], order=3)],
				["S_pd_meas_err", [], "continuous", .01],  #mA/W
				#snr
				["tau", [], "continuous", 1800],
				["skyspec", [], "object", np.genfromtxt(skyfile_path,usecols=[0,1],names=['waves_nm','skyflux'])],
				#llamas
				["resolution", [], "continuous", 2200.0],
				["wave_min", [], "continuous", _wave_min],
				["wave_max", [], "continuous", _wave_max],
				["wave_redgreen", [], "continuous", 690.0],
				["wave_greenblue", [], "continuous", 480.0],
				["f_collimator", [], "continuous", 200.0],
				["f_camera", [], "continuous", 70.0],
				["fratio_collimator", [], "continuous", 4.0],
				["fratio_camera", [], "continuous", 1.28],
				["blade_obscure", [], "continuous", 0.9],
				["d_beam", [], "continuous", 50.0],
				#fiber
				["fiber_theta", [], "continuous", 0.75],
				["fiber_dcore", [], "continuous", 110],
				["microlens", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
				["fiber_ar", [], "functional", define_functional_from_file(_dir+"fiber_ar.txt")],
				["fiber_internal", [], "functional", define_functional_from_file(_dir+"Polymicro_FBPI_8m.txt")],
				["frd_meas_err", [], "continuous", 0.068], #analysis based on test measuring a few fibers
				#spectrograph
				["spect_model_tracker", [], "object", spectrograph_models],
				["collimator", [], "functional", define_functional_from_file(_dir+"dielectric_mirror.txt")],
				["prism", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
				["lens1", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
				["lens2", [], "functional", define_functional_from_file(_dir+"ECI_PBM8Y.txt")],
				["lens3", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
				["lens4", [], "functional", define_functional_from_file(_dir+"ECI_PBM8Y.txt")],
				["lens5", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
				["lens6", [], "functional", define_functional_from_file(_dir+"ECI_PBM8Y.txt")],
				["lens7", [], "functional", define_functional_from_file(_dir+"ECI_PBM8Y.txt")],
				["sensor_window", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
				["sensor_glass_red", [], "functional", define_functional_from_file(_dir+"llamas_internal_red.txt")],
				["sensor_glass_gre", [], "functional", define_functional_from_file(_dir+"llamas_internal.txt")],
				["sensor_glass_blu", [], "functional", define_functional_from_file(_dir+"llamas_internal_blue.txt")],
				#optical measurements
				["vph_meas_stddev", [], "continuous", .001],
				["sl_meas_stddev", [], "continuous", .001],
				["bg_meas_stddev", [], "continuous", .001],
				#cost
				["testbed_setup", ["nonrandom", [1800]], "continuous", 1800], #WAG
				#["C_engineer", ["nonrandom", [0.00694444444]], "continuous", 0.00694444444] #WAG $/s, from $25/hr
				["C_engineer", ["nonrandom", [1]], "continuous", 1] #just count time
			]


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
eta = snr_likelihood_fn
H = sensitivity_hlva
Gamma = cost_model
llamas_snr = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs, x_defs)


def update_llamas_problem(llamas_snr, d):
	d_masked = [(math.floor(dd) if llamas_snr.d_masks[i]=='discrete' else dd) for i,dd in enumerate(d)]
	d_dict = dict(zip(llamas_snr.d_names, d_masked))
	num_vph = d_dict["d_vph_n_pts"]
	num_dichroic = d_dict["d_dichroic_n_pts"]
	
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
		eta, 
		H, 
		Gamma, 
		theta_defs, 
		mod_y_defs, #!
		d_defs, 
		x_defs
	)
		
	return llamas_snr_d