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
from approx.regression_models import *
#llamas
from problems.llamas_snr.snr_exp_models import *
from problems.llamas_snr.snr_system_model import *

#set path variables
_dir = "./llamas-etc/COATINGS/"

_wave_min = 350.0
_wave_max = 975.0
_wave_redgreen = 690.0
_wave_greenblue = 480.0
_bandpass = _wave_max - _wave_min
#my physical argument here is that we should consider measurements at these distances apart to be essentially independent
#and this is based on the approximate spectral resolution of the instrument
#R = lambda / dlambda = 2200
#so for lambda=350, dlambda=0.159
#and for lambda=975, dlambda=0.443
_lengthscale = 1.0
thru_param_var = .01**2

prior_gain_SN1 = ["gamma_mv",  [0.999,0.2**2]] #mean, variance
prior_rn_SN1 = ["gamma_mv", [2.32,0.25**2]]
prior_dc_SN1 = ["gamma_mv", [0.00238,.001**2]]

prior_gain_SN3 = ["gamma_mv",  [1.008,0.2**2]] #mean, variance
prior_rn_SN3 = ["gamma_mv", [2.35,0.25**2]]
prior_dc_SN3 = ["gamma_mv", [0.00267,.001**2]]

#ppts, meanfn = get_ppts_meanfn_file(_dir+"CCD42-40_dd.txt", 3)
#prior_qe_red = ["gp_expquad", [.05, _lengthscale, ppts, meanfn]]
#ppts, meanfn = get_ppts_meanfn_file(_dir+"CCD42-40_green.txt", 3)
#prior_qe_gre = ["gp_expquad", [.05, _lengthscale, ppts, meanfn]]
#ppts, meanfn = get_ppts_meanfn_file(_dir+"CCD42-40_blue.txt", 3)
#prior_qe_blu = ["gp_expquad", [.05, _lengthscale, ppts, meanfn]]

coeffs, inter, _, _ = linreg_fourier_throughput_file(_dir+"CCD42-40_dd.txt", 2, doPlot=False, doErr=True)
prior_qe_red_t0 = ["gaussian", [inter,thru_param_var]] 
prior_qe_red_t1 = ["gaussian", [coeffs[0],thru_param_var]]
prior_qe_red_t2 = ["gaussian", [coeffs[1],thru_param_var]]
prior_qe_red_t3 = ["gaussian", [coeffs[2],thru_param_var]]
prior_qe_red_t4 = ["gaussian", [coeffs[3],thru_param_var]]
coeffs, inter, _, _ = linreg_fourier_throughput_file(_dir+"CCD42-40_green.txt", 2, doPlot=False, doErr=True)
prior_qe_gre_t0 = ["gaussian", [inter,thru_param_var]]
prior_qe_gre_t1 = ["gaussian", [coeffs[0],thru_param_var]]
prior_qe_gre_t2 = ["gaussian", [coeffs[1],thru_param_var]]
prior_qe_gre_t3 = ["gaussian", [coeffs[2],thru_param_var]]
prior_qe_gre_t4 = ["gaussian", [coeffs[3],thru_param_var]]
coeffs, inter, _, _ = linreg_fourier_throughput_file(_dir+"CCD42-40_blue.txt", 2, doPlot=False, doErr=True)
prior_qe_blu_t0 = ["gaussian", [inter,thru_param_var]]
prior_qe_blu_t1 = ["gaussian", [coeffs[0],thru_param_var]]
prior_qe_blu_t2 = ["gaussian", [coeffs[1],thru_param_var]]
prior_qe_blu_t3 = ["gaussian", [coeffs[2],thru_param_var]]
prior_qe_blu_t4 = ["gaussian", [coeffs[3],thru_param_var]]

#ppts, meanfn = get_ppts_meanfn_file(_dir+"wasach_llamas2200_red.txt", 2)
#prior_gp_vph_red = ["gp_expquad", [.0267, _lengthscale, ppts, meanfn]]
#ppts, meanfn = get_ppts_meanfn_file(_dir+"wasach_llamas2200_green.txt", 2)
#prior_gp_vph_gre = ["gp_expquad", [.0267, _lengthscale, ppts, meanfn]]
#ppts, meanfn = get_ppts_meanfn_file(_dir+"wasach_llamas2200_blue.txt", 2)
#prior_gp_vph_blu = ["gp_expquad", [.0267, _lengthscale, ppts, meanfn]]

a, b, c, d = p3_fit_throughput_file(_dir+"wasach_llamas2200_red.txt", doPlot=False, doErr=True)
prior_vph_red_t0 = ["gaussian", [a,thru_param_var]] 
prior_vph_red_t1 = ["gaussian", [b,thru_param_var]]
prior_vph_red_t2 = ["gaussian", [c,thru_param_var]]
prior_vph_red_t3 = ["gaussian", [d,thru_param_var]]
a, b, c, d = p3_fit_throughput_file(_dir+"wasach_llamas2200_green.txt", doPlot=False, doErr=True)
prior_vph_gre_t0 = ["gaussian", [a,thru_param_var]]
prior_vph_gre_t1 = ["gaussian", [b,thru_param_var]]
prior_vph_gre_t2 = ["gaussian", [c,thru_param_var]]
prior_vph_gre_t3 = ["gaussian", [d,thru_param_var]]
a, b, c, d = p3_fit_throughput_file(_dir+"wasach_llamas2200_blue.txt", doPlot=False, doErr=True)
prior_vph_blu_t0 = ["gaussian", [a,thru_param_var]]
prior_vph_blu_t1 = ["gaussian", [b,thru_param_var]]
prior_vph_blu_t2 = ["gaussian", [c,thru_param_var]]
prior_vph_blu_t3 = ["gaussian", [d,thru_param_var]]

#ppts, meanfn = get_ppts_meanfn_file(_dir+"ECI_FusedSilica_sl_prior.txt", 3)
#prior_gp_sl = ["gp_expquad", [.1, _lengthscale, ppts, meanfn]]
#ppts, meanfn = get_ppts_meanfn_file(_dir+"ECI_FusedSilica_bg_prior.txt", 3)
#prior_gp_bg = ["gp_expquad", [.1, _lengthscale, ppts, meanfn]]

lval, steppt, rval, power = sigmoid_fit_throughput_file(_dir+"ECI_FusedSilica_sl_prior.txt", doPlot=False, doErr=True)
prior_sl_t0 = ["gaussian", [lval,thru_param_var]]
prior_sl_t1 = ["gaussian", [steppt,thru_param_var]]
prior_sl_t2 = ["gaussian", [rval,thru_param_var]]
prior_sl_t3 = ["gaussian", [power,thru_param_var]]
lval, steppt, rval, power = sigmoid_fit_throughput_file(_dir+"ECI_FusedSilica_bg_prior.txt", doPlot=False, doErr=True)
prior_bg_t0 = ["gaussian", [lval,thru_param_var]]
prior_bg_t1 = ["gaussian", [steppt,thru_param_var]]
prior_bg_t2 = ["gaussian", [rval,thru_param_var]]
prior_bg_t3 = ["gaussian", [power,thru_param_var]]

prior_frd = ["gamma_mv", [0.077,0.022**2]]


#these priors are based on requirements that were met, see Camera Qual Report
theta_defs = [                             #mean, variance
					["gain_red", prior_gain_SN1, "continuous"],
					["gain_gre", prior_gain_SN3, "continuous"],
					["gain_blu", prior_gain_SN3, "continuous"],
					["rn_red", prior_rn_SN1, "continuous"],
					["rn_gre", prior_rn_SN3, "continuous"],
					["rn_blu", prior_rn_SN3, "continuous"],
					["dc_red", prior_dc_SN1, "continuous"],
					["dc_gre", prior_dc_SN3, "continuous"],
					["dc_blu", prior_dc_SN3, "continuous"],
					#["qe_red", prior_qe_red, "functional"],
					["qe_red_t0", prior_qe_red_t0, "continuous"],
					["qe_red_t1", prior_qe_red_t1, "continuous"],
					["qe_red_t2", prior_qe_red_t2, "continuous"],
					["qe_red_t3", prior_qe_red_t3, "continuous"],
					["qe_red_t4", prior_qe_red_t4, "continuous"],
					#["qe_gre", prior_qe_gre, "functional"],
					["qe_gre_t0", prior_qe_gre_t0, "continuous"],
					["qe_gre_t1", prior_qe_gre_t1, "continuous"],
					["qe_gre_t2", prior_qe_gre_t2, "continuous"],
					["qe_gre_t3", prior_qe_gre_t3, "continuous"],
					["qe_gre_t4", prior_qe_gre_t4, "continuous"],
					#["qe_blu", prior_qe_blu, "functional"],
					["qe_blu_t0", prior_qe_blu_t0, "continuous"],
					["qe_blu_t1", prior_qe_blu_t1, "continuous"],
					["qe_blu_t2", prior_qe_blu_t2, "continuous"],
					["qe_blu_t3", prior_qe_blu_t3, "continuous"],
					["qe_blu_t4", prior_qe_blu_t4, "continuous"],
					#["vph_thru_red", prior_gp_vph_red, "functional"],
					["vph_red_t0", prior_vph_red_t0, "continuous"],
					["vph_red_t1", prior_vph_red_t1, "continuous"],
					["vph_red_t2", prior_vph_red_t2, "continuous"],
					["vph_red_t3", prior_vph_red_t3, "continuous"],
					#["vph_thru_gre", prior_gp_vph_gre, "functional"],
					["vph_gre_t0", prior_vph_gre_t0, "continuous"],
					["vph_gre_t1", prior_vph_gre_t1, "continuous"],
					["vph_gre_t2", prior_vph_gre_t2, "continuous"],
					["vph_gre_t3", prior_vph_gre_t3, "continuous"],
					#["vph_thru_blu", prior_gp_vph_blu, "functional"],
					["vph_blu_t0", prior_vph_blu_t0, "continuous"],
					["vph_blu_t1", prior_vph_blu_t1, "continuous"],
					["vph_blu_t2", prior_vph_blu_t2, "continuous"],
					["vph_blu_t3", prior_vph_blu_t3, "continuous"],
					#["sl_thru_dichroic", prior_gp_sl, "functional"],
					["sl_t0", prior_sl_t0, "continuous"],
					["sl_t1", prior_sl_t1, "continuous"],
					["sl_t2", prior_sl_t2, "continuous"],
					["sl_t3", prior_sl_t3, "continuous"],
					#["bg_thru_dichroic", prior_gp_bg, "functional"],
					["bg_t0", prior_bg_t0, "continuous"],
					["bg_t1", prior_bg_t1, "continuous"],
					["bg_t2", prior_bg_t2, "continuous"],
					["bg_t3", prior_bg_t3, "continuous"],
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
				#"y_qe_red_pts",
				"y_qe_red_t0", 
				"y_qe_red_t1", 
				"y_qe_red_t2", 
				"y_qe_red_t3", 
				"y_qe_red_t4", 
				#"y_qe_gre_pts",
				"y_qe_gre_t0", 
				"y_qe_gre_t1", 
				"y_qe_gre_t2", 
				"y_qe_gre_t3", 
				"y_qe_gre_t4", 
				#"y_qe_blu_pts",
				"y_qe_blu_t0", 
				"y_qe_blu_t1", 
				"y_qe_blu_t2", 
				"y_qe_blu_t3", 
				"y_qe_blu_t4", 
				#"y_vph_red_pts",
				"y_vph_red_t0", 
				"y_vph_red_t1", 
				"y_vph_red_t2", 
				"y_vph_red_t3", 
				#"y_vph_gre_pts",
				"y_vph_gre_t0", 
				"y_vph_gre_t1", 
				"y_vph_gre_t2", 
				"y_vph_gre_t3", 
				#"y_vph_blu_pts",
				"y_vph_blu_t0", 
				"y_vph_blu_t1", 
				"y_vph_blu_t2", 
				"y_vph_blu_t3", 
				#"y_sl_pts",
				"y_sl_t0", 
				"y_sl_t1", 
				"y_sl_t2", 
				"y_sl_t3", 
				#"y_bg_pts",
				"y_bg_t0", 
				"y_bg_t1", 
				"y_bg_t2", 
				"y_bg_t3", 
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
with fits.open('./llamas-etc/'+'lbg_shapley.fits') as hdul:
	shapley_wv   = hdul[2].data * (1+2.5) / 10.0
	shapley_spec = hdul[0].data

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
				["spectral_power", [], "continuous", 1e-4], #watts / nm, Energetiq LDLS, Krishnamurthy et al. 2016
				#snr
				["tau", [], "continuous", 5*3600],
				["skyspec", [], "object", np.genfromtxt(skyfile_path,usecols=[0,1],names=['waves_nm','skyflux'])],
				["shapley_wv", [], "object", shapley_wv],
				["shapley_spec", [], "object", shapley_spec],
				#llamas
				["resolution", [], "continuous", 2200.0],
				["wave_min", [], "continuous", _wave_min],
				["wave_max", [], "continuous", _wave_max],
				["wave_redgreen", [], "continuous", _wave_redgreen],
				["wave_greenblue", [], "continuous", _wave_greenblue],
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