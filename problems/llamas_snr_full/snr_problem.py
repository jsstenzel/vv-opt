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
from approx.gaussian_process import *
from approx.learn_gp import *
from approx.regression_models import *
#llamas
from problems.llamas_snr_full.snr_exp_models import *
from problems.llamas_snr_full.snr_system_model import *
from problems.llamas_snr_full.snr_cost_model import *

def learn_llamas_snr_priors(verbose=True, doPlot=True, doSave=True):
	#Directories to find prior training data
	_dir = "./llamas-etc/COATINGS/"
	_reddir = "./vendor data/lenses/camera_RED/"
	_gredir = "./vendor data/lenses/camera_GREEN/"
	_bludir = "./vendor data/lenses/camera_BLUE/"
	_qedir = "./vendor data/quantum efficiency/"
	
	#Parameters for prior training
	#T_std = 1 #This made really nice plots, but isnt physically motivated
	T_std = 0.0001 #is more physically correct, dT = 0.01%
	QE_std = 0.01
	
	#Directory to save new prior def files to
	sloc = './priors/'

	###Quantum efficiency priors
	learn_save_gp_prior_to_file(sloc+'prior_gp_qe_red', [_qedir+"qe_red_basic_nir.txt"], QE_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gp_qe_gre', [_qedir+"qe_green_basic_midband.txt"], QE_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gp_qe_blu', [_qedir+"qe_blue_enhanced_broadband.txt"], QE_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)

	###VPH priors
	learn_save_gp_prior_to_file(sloc+'prior_gp_vph_red', [_dir+"wasach_llamas2200_red.txt"], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gp_vph_gre', [_dir+"wasach_llamas2200_green.txt"], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gp_vph_blu', [_dir+"wasach_llamas2200_blue.txt"], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	
	###Red Lens priors -- manufacturer curves for 60-30110
	learn_save_gp_prior_to_file(sloc+'prior_red1', [_reddir+"red1_2241C.txt", _reddir+'red1_2243C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_red2', [_reddir+"red2_2209C.txt", _reddir+'red2_2212C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_red3', [_reddir+"red3_2205C.txt", _reddir+'red3_2208C.txt', _reddir+'red3_2209C.txt', _reddir+'red3_2212C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_red4', [_reddir+"red4_2209C.txt", _reddir+'red4_2214C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_red5', [_reddir+'red5_2217C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_red6', [_reddir+'red6_9671C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_red7', [_reddir+'red7_2241C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	
	###Green Lens priors -- manufacturer curves for 60-30110
	learn_save_gp_prior_to_file(sloc+'prior_gre1', [_gredir+"green1_8992C.txt", _gredir+'green1_8995C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gre2', [_gredir+"green2_9455C.txt", _gredir+'green2_9459C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gre3', [_gredir+"green3_9463C.txt", _gredir+'green3_9466C.txt', _gredir+'green3_9468C.txt', _gredir+'green3_9480C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gre4', [_gredir+"green4_9455C.txt", _gredir+'green4_9459C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gre5', [_gredir+'green5_9476C.txt'], T_std, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gre6', [_gredir+'green6_8970C.txt', _gredir+'green6_9506C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_gre7', [_gredir+'green7_9478C.txt'], T_std, doPlot=doPlot, doPrint=verbose, careful=True)
	
	###Blue Lens priors -- manufacturer curves for 60-30110
	learn_save_gp_prior_to_file(sloc+'prior_blu1', [_bludir+"blue1_9371C.txt", _bludir+'blue1_9376C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_blu2', [_bludir+"blue2_1961C.txt", _bludir+'blue2_1965C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_blu3', [_bludir+"blue3_9360C.txt", _bludir+'blue3_9363C.txt', _bludir+'blue3_9371C.txt', _bludir+'blue3_9366C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)

	learn_save_gp_prior_to_file(sloc+'prior_blu4', [_bludir+"blue4_9363C.txt", _bludir+'blue4_9373C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_blu5', [_bludir+'blue5_1961C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_blu6', [_bludir+'blue6_1961C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_blu7', [_bludir+'blue7_9360C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	learn_save_gp_prior_to_file(sloc+'prior_blu8', [_bludir+'blue8_9412C.txt', _bludir+'blue8_9414C.txt'], T_std, save=doSave, doPlot=doPlot, doPrint=verbose, careful=True)
	

def construct_llamas_snr_problem(verbose_probdef=False):
	print("Constructing llamas_snr_full priors ...",flush=True)

	#set path variables
	_dir = "./llamas-etc/COATINGS/"
	
	#Directory to load prior def files from
	sloc = './priors/'

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
	_lengthscale = 50.0

	###Camera priors
	prior_gain_SN1 = ["gamma_mv",  [0.999,0.2**2]] #mean, variance
	prior_rn_SN1 = ["gamma_mv", [2.32,0.25**2]]
	prior_dc_SN1 = ["gamma_mv", [0.00238,.001**2]]

	prior_gain_SN3 = ["gamma_mv",  [1.008,0.2**2]] #mean, variance
	prior_rn_SN3 = ["gamma_mv", [2.35,0.25**2]]
	prior_dc_SN3 = ["gamma_mv", [0.00267,.001**2]]

	###Quantum efficiency priors
	prior_gp_qe_red = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gp_qe_red',doPlot=verbose_probdef))]
	prior_gp_qe_gre = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gp_qe_gre',doPlot=verbose_probdef))]
	prior_gp_qe_blu = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gp_qe_blu',doPlot=verbose_probdef))]

	###VPH priors
	prior_gp_vph_red = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gp_vph_red',doPlot=verbose_probdef))]
	prior_gp_vph_gre = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gp_vph_gre',doPlot=verbose_probdef))]
	prior_gp_vph_blu = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gp_vph_blu',doPlot=verbose_probdef))]
	
	###Dichroic priors
	#TODO gp_expquad is not really a suitable kernel for learning a GP prior, due to the big discontinuity
	#So for now, I will continue to define the mean function and kernel properties manually
	ppts, _ = get_ppts_meanfn_file(_dir+"ECI_FusedSilica_sl_prior.txt", 3, doPlot=False)
	lval_sl, steppt_sl, rval_sl, power_sl, _ = sigmoid_fit_throughput_file(_dir+"ECI_FusedSilica_sl_prior.txt", doPlot=verbose_probdef, doErr=False)
	meanfn_sl_prior = t_to_u(throughput_from_sigmoidfit_coeffs(lval_sl, steppt_sl, rval_sl, power_sl, ppts))
	var_u = 1.0 #kind of a WAG
	prior_gp_sl = ["gp_expquad", [var_u, _lengthscale, ppts, meanfn_sl_prior]]

	ppts, _ = get_ppts_meanfn_file(_dir+"ECI_FusedSilica_bg_prior.txt", 3, doPlot=False)
	lval_bg, steppt_bg, rval_bg, power_bg, _ = sigmoid_fit_throughput_file(_dir+"ECI_FusedSilica_bg_prior.txt", doPlot=verbose_probdef, doErr=False)
	meanfn_bg_prior = t_to_u(throughput_from_sigmoidfit_coeffs(lval_bg, steppt_bg, rval_bg, power_bg, ppts))
	var_u = 1.0 #kind of a WAG
	prior_gp_bg = ["gp_expquad", [var_u, _lengthscale, ppts, meanfn_bg_prior]]

	###Collimator priors, TODO I can do better later
	coll_prior_pts = list(np.arange(325,975,1))
	coll_mean_pts = [0.99 for _ in coll_prior_pts]
	coll_prior_pts.append(1000)
	coll_mean_pts.append(0.98)
	var_u = 1.0 #kind of a WAG
	ls = 8
	coll_fn = t_to_u(coll_mean_pts)
	prior_gp_coll = ["gp_expquad", [var_u, ls, coll_prior_pts, coll_fn]]
	
	###Red Lens priors -- manufacturer curves for 60-30110
	prior_red1 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_red1',doPlot=verbose_probdef))]	
	prior_red2 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_red2',doPlot=verbose_probdef))]
	prior_red3 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_red3',doPlot=verbose_probdef))]	
	prior_red4 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_red4',doPlot=verbose_probdef))]	
	prior_red5 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_red5',doPlot=verbose_probdef))]	
	prior_red6 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_red6',doPlot=verbose_probdef))]	
	prior_red7 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_red7',doPlot=verbose_probdef))]
	
	###Green Lens priors -- manufacturer curves for 60-30110
	prior_gre1 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gre1',doPlot=verbose_probdef))]	
	prior_gre2 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gre2',doPlot=verbose_probdef))]
	prior_gre3 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gre3',doPlot=verbose_probdef))]	
	prior_gre4 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gre4',doPlot=verbose_probdef))]	
	prior_gre5 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gre5',doPlot=verbose_probdef))]	
	prior_gre6 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gre6',doPlot=verbose_probdef))]	
	prior_gre7 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_gre7',doPlot=verbose_probdef))]
	
	###Blue Lens priors -- manufacturer curves for 60-30110
	prior_blu1 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu1',doPlot=verbose_probdef))]
	prior_blu2 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu2',doPlot=verbose_probdef))]
	prior_blu3 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu3',doPlot=verbose_probdef))]
	prior_blu4 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu4',doPlot=verbose_probdef))]
	prior_blu5 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu5',doPlot=verbose_probdef))]
	prior_blu6 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu6',doPlot=verbose_probdef))]
	prior_blu7 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu7',doPlot=verbose_probdef))]
	prior_blu8 = ["gp_expquad", list(learn_load_gp_prior_from_file(sloc+'prior_blu8',doPlot=verbose_probdef))]

	###Fiber priors
	prior_frd = ["beta", [3.434339623321249, 133.9392453095287]]

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
						["qe_red_t", prior_gp_qe_red, "continuous"],
						["qe_gre_t", prior_gp_qe_gre, "continuous"],
						["qe_blu_t", prior_gp_qe_blu, "continuous"],
						["vph_red_t", prior_gp_vph_red, "continuous"],
						["vph_gre_t", prior_gp_vph_gre, "continuous"],
						["vph_blu_t", prior_gp_vph_blu, "continuous"],
						["sl_t", prior_gp_sl, "continuous"],
						["bg_t", prior_gp_bg, "continuous"],
						["coll_t", prior_gp_coll, "continuous"],
						["red_l1_t", prior_red1, "continuous"],
						["red_l2_t", prior_red2, "continuous"],
						["red_l3_t", prior_red3, "continuous"],
						["red_l4_t", prior_red4, "continuous"],
						["red_l5_t", prior_red5, "continuous"],
						["red_l6_t", prior_red6, "continuous"],
						["red_l7_t", prior_red7, "continuous"],
						["gre_l1_t", prior_gre1, "continuous"],
						["gre_l2_t", prior_gre2, "continuous"],
						["gre_l3_t", prior_gre3, "continuous"],
						["gre_l4_t", prior_gre4, "continuous"],
						["gre_l5_t", prior_gre5, "continuous"],
						["gre_l6_t", prior_gre6, "continuous"],
						["gre_l7_t", prior_gre7, "continuous"],
						["blu_l1_t", prior_blu1, "continuous"],
						["blu_l2_t", prior_blu2, "continuous"],
						["blu_l3_t", prior_blu3, "continuous"],
						["blu_l4_t", prior_blu4, "continuous"],
						["blu_l5_t", prior_blu5, "continuous"],
						["blu_l6_t", prior_blu6, "continuous"],
						["blu_l7_t", prior_blu7, "continuous"],
						["blu_l8_t", prior_blu8, "continuous"],
						["fiber_frd", prior_frd, "continuous"]
					]


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
					"y_qe_red_t0", 
					"y_qe_red_t1", 
					"y_qe_red_t2", 
					"y_qe_red_t3", 
					"y_qe_red_t4", 
					"y_qe_gre_t0", 
					"y_qe_gre_t1", 
					"y_qe_gre_t2", 
					"y_qe_gre_t3", 
					"y_qe_gre_t4", 
					"y_qe_blu_t0", 
					"y_qe_blu_t1", 
					"y_qe_blu_t2", 
					"y_qe_blu_t3", 
					"y_qe_blu_t4", 
					"y_vph_red_p0", 
					"y_vph_red_p1", 
					"y_vph_red_p2", 
					"y_vph_gre_p0", 
					"y_vph_gre_p1", 
					"y_vph_gre_p2",
					"y_vph_blu_p0", 
					"y_vph_blu_p1", 
					"y_vph_blu_p2", 
					"y_sl_t0", 
					"y_sl_t1", 
					"y_sl_t2", 
					"y_sl_t3", 
					"y_bg_t0", 
					"y_bg_t1", 
					"y_bg_t2", 
					"y_bg_t3", 
					"y_coll_t",
					"y_red_cam_t",
					"y_gre_cam_t",
					"y_blu_cam_t",
					"y_frd"
				]


	d_defs_full = [
					["t_gain", ['uniform', [.1, 600]], "continuous"], #gain
					["I_gain", ['uniform', [0, 100]], "discrete"],    #gain
					["n_meas_rn", ['uniform', [0, 50]], "discrete"],  #rn
					["d_num", ['uniform', [2, 25]], "discrete"],      #dc
					["d_max", ['uniform', [0, 12000]], "continuous"], #dc
					["d_pow", ['uniform', [0.01,3]], "continuous"],      #dc
					
					["n_qe", ['uniform', [0, 100]], "discrete"],   #qe
					["t_qe", ['uniform', [.1, 600]], "continuous"],#qe
					
					["d_vph_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_dichroic_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_coll_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_redcam_n_pts", ['uniform', [0,(_wave_max-_wave_redgreen)*10]], "discrete"], #i guess??
					["d_greencam_n_pts", ['uniform', [0,(_wave_redgreen-_wave_greenblue)*10]], "discrete"], #i guess??
					["d_bluecam_n_pts", ['uniform', [0,(_wave_greenblue-_wave_min)*10]], "discrete"], #i guess??
					["d_frd_n_meas", ['uniform', [0,2400]], "discrete"],
				]
	d_defs_non0 = [
					["t_gain", ['uniform', [.1, 600]], "continuous"], #gain
					["I_gain", ['uniform', [1, 100]], "discrete"],    #gain
					["n_meas_rn", ['uniform', [1, 50]], "discrete"],  #rn
					["d_num", ['uniform', [2, 25]], "discrete"],      #dc
					["d_max", ['uniform', [1, 12000]], "continuous"], #dc
					["d_pow", ['uniform', [0.01,3]], "continuous"],      #dc
					
					["n_qe", ['uniform', [1, 100]], "discrete"],   #qe
					["t_qe", ['uniform', [.1, 600]], "continuous"],#qe
					
					["d_vph_n_pts", ['uniform', [1,_bandpass*10]], "discrete"],
					["d_dichroic_n_pts", ['uniform', [1,_bandpass*10]], "discrete"],
					["d_coll_n_pts", ['uniform', [1,_bandpass*10]], "discrete"],
					["d_redcam_n_pts", ['uniform', [1,(_wave_max-_wave_redgreen)*10]], "discrete"], #i guess??
					["d_greencam_n_pts", ['uniform', [1,(_wave_redgreen-_wave_greenblue)*10]], "discrete"], #i guess??
					["d_bluecam_n_pts", ['uniform', [1,(_wave_greenblue-_wave_min)*10]], "discrete"], #i guess??
					["d_frd_n_meas", ['uniform', [1,2400]], "discrete"],
				]
	d_defs_discrete = [
					["t_gain", ['uniform', [0, 600]], "discrete"], #gain
					["I_gain", ['uniform', [0, 100]], "discrete"],    #gain
					["n_meas_rn", ['uniform', [0, 50]], "discrete"],  #rn
					["d_num", ['uniform', [2, 25]], "discrete"],      #dc
					["d_max", ['uniform', [.1, 12000]], "discretized100"], #dc
					["d_pow", ['uniform', [0.01,3]], "discretized100"],      #dc
					
					["n_qe", ['uniform', [0, 100]], "discrete"],   #qe
					["t_qe", ['uniform', [0, 600]], "discrete"],#qe
					
					["d_vph_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_dichroic_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_coll_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_redcam_n_pts", ['uniform', [0,(_wave_max-_wave_redgreen)*10]], "discrete"], #i guess??
					["d_greencam_n_pts", ['uniform', [0,(_wave_redgreen-_wave_greenblue)*10]], "discrete"], #i guess??
					["d_bluecam_n_pts", ['uniform', [0,(_wave_greenblue-_wave_min)*10]], "discrete"], #i guess??
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
					#rn exp
					["t_rn", [], "continuous", .1], #100 ms exposure
					#dc exp
					["t_0", [], "continuous", 0.1], #100ms baseline exposure assumed
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
					["prism", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
					["sensor_window", [], "functional", define_functional_from_file(_dir+"ECI_FusedSilica.txt")],
					["sensor_glass_red", [], "functional", define_functional_from_file(_dir+"llamas_internal_red.txt", order=2)],
					["sensor_glass_gre", [], "functional", define_functional_from_file(_dir+"llamas_internal.txt", order=2)],
					["sensor_glass_blu", [], "functional", define_functional_from_file(_dir+"llamas_internal_blue.txt", order=2)],
					#optical measurements
					["vph_meas_stddev", [], "continuous", .002], #Measurement Considerations When Specifying Optical Coatings, Pete Kupinski and Angus Macleod. This paper indicates a best case +- 0.1% T for commercial measurements of highly transmissive coatings. Lets add margin by doubling that.
					["sl_meas_stddev", [], "continuous", .002], #ibid.
					["bg_meas_stddev", [], "continuous", .002], #ibid.
					["coll_meas_stddev", [], "continuous", .002], #ibid.
					["lens_meas_err", [], "continuous", .002], #ibid.
					#cost
					["t_gain_setup", [], "continuous", 1200], #rough estimate based on experience
					["t_gain_buffer", [], "continuous", 5], #rough estimate based on experience
					["t_rn_buffer", [], "continuous", 5], #rough estimate based on experience
					["t_dc_buffer", [], "continuous", 5], #rough estimate based on experience
					["fp_testbed_setup", [], "continuous", 1800], #rough estimate based on experience
					["t_qe_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
					["t_qe_buffer", [], "continuous", 60], #WAG - need to adjust photodiode or sensor or something, so say about a minute
					["t_vph_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
					["t_vph_per_pt", [], "continuous", 5*60], #WAG - individual measurements -5 min each with buffer
					["t_dichroic_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
					["t_dichroic_per_pt", [], "continuous", 0.01], #WAG - this captures tradeoff between spectral resolution and integration time
					["t_coll_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
					["t_coll_per_pt", [], "continuous", 0.05], #WAG - this captures tradeoff between spectral resolution and integration time
					["t_camera_test_setup", [], "continuous", 1800], #WAG this setup is probably more complicated
					["t_camera_per_pt", [], "continuous", 5*60], #WAG - individual measurements - 5 min each with buffer
					["t_frd_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
					["t_frd_test", [], "continuous", 600], #WAG i think this test was probably pretty fiddly, since we only tested 10
					#["C_engineer", [], "continuous", 0.00694444444] #WAG $/s, from $25/hr
					["C_engineer", [], "continuous", 1] #just count time
				]


	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	eta = snr_likelihood_fn
	H = sensitivity_hlva
	Gamma = snr_cost
	llamas_snr = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs_discrete, x_defs)
	print("llamas_snr_full problem constructed.",flush=True)
	return llamas_snr
	
if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--train_priors', action='store_true', required=False, help='Run the learn_llamas_snr_priors function')
	args = parser.parse_args()
	
	if args.train_priors:
		learn_llamas_snr_priors(verbose=True)
	else:
		llamas_snr = construct_llamas_snr_problem(verbose_probdef=True)
		print(llamas_snr)