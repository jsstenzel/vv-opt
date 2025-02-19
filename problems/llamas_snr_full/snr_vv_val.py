import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv
import dill

sys.path.append('../..')
#llamas
from problems.llamas_snr_full.snr_problem import *
#analysis
from obed.obed_multivar import *
from obed.obed_gbi import *
#from obed.pdf_estimation import *
from inference.bn_modeling import *
from uq.uncertainty_propagation import *
from uq.sensitivity_analysis import *
from uq.saltelli_gsa import *
from uq.gsa_convergence import *
from opt.ngsa import *

################################
#Useful definitions
################################

################################
#Analysis functions
################################

def get_y_hist(problem):
	_dir = "./llamas-etc/COATINGS/"
	_bgdir = "./historical data/BG dichroic/"
	_sldir = "./historical data/SL dichroic/"
	_qedir = "./vendor data/quantum efficiency/"	
	
	###get X, other definitions
	x = problem.x_dict
	prior_mean = problem.prior_mean_dict 
	_wave_min = 350.0
	_wave_max = 975.0
	_wave_redgreen = 690.0
	_wave_greenblue = 480.0
	_bandpass = _wave_max - _wave_min
	red_max = _wave_max
	red_min = _wave_redgreen
	gre_max = _wave_redgreen
	gre_min = _wave_greenblue
	blu_max = _wave_greenblue
	blu_min = _wave_min
	
	###CCD
	#["SN20012", "Red 3", "nodither", [2.44,.025], [0.00113,.000248], [1.031,.0106]],
	#["SN20014", "Green 4", "dither", [2.54,0.0262], [0.00147,0.000131], [1.061,0.0109],
	#["SN20019", "Blue 6", "dither", [2.26,.023], [0.00054,2.34e-05], [1.025,.0104]],
	gain_red = 1.031
	gain_gre = 1.061
	gain_blu = 1.025
	y_rn_red = 2.44
	y_rn_gre = 2.54
	y_rn_blu = 2.26
	y_dc_red = 0.00113
	y_dc_gre = 0.00147
	y_dc_blu = 0.00054
	#y_gain isnt exactly the same dimension at theta_gain, unlike the other two experiments, so we gotta convert it
	y_gain_red = gain_exp(gain_red, 0, 0, 20, 30, x, prior_mean, err=False)
	y_gain_gre = gain_exp(gain_gre, 0, 0, 20, 30, x, prior_mean, err=False)
	y_gain_blu = gain_exp(gain_blu, 0, 0, 20, 30, x, prior_mean, err=False)
	
	###for QE, since historical design didnt measure it, 
	#we have to simulate the experiment under d=0 condition with the prior mean. Awkward imputation but ok
	y_qe_red_t = quantum_efficiency_exp(None, y_gain_red, y_rn_red, 0, 0, red_min, red_max, red_max-blu_min, x, prior_mean["qe_red_t"], False)
	y_qe_gre_t = quantum_efficiency_exp(None, y_gain_gre, y_rn_gre, 0, 0, gre_min, gre_max, red_max-blu_min, x, prior_mean["qe_gre_t"], False)
	y_qe_blu_t = quantum_efficiency_exp(None, y_gain_blu, y_rn_blu, 0, 0, blu_min, blu_max, red_max-blu_min, x, prior_mean["qe_blu_t"], False)
	
	###VPH - idk which are in which spectrograph, just picking out of a hat
	"""
	6414G-07-02	475nm	55.71%
	6414G-07-02	550nm	75.20%
	6414G-07-02	675nm	56.60%

	6414R-04-01	700nm	85.68%
	6414R-04-01	825nm	84.15%
	6414R-04-01	925nm	65.05%

	6414B-09-01	375nm	66.69%
	6414B-09-01	425nm	76.17%
	6414B-09-01	475nm	57.58%
	"""
	#make the thru curves from here
	y_vph_red_p, _ = poly_fit_throughput([700, 825, 925], [0.8568, 0.8415, 0.6505], 2, doPlot=False, doErr=False, doCov=False)
	y_vph_gre_p, _ = poly_fit_throughput([475, 550, 675], [0.5571, 0.7520, 0.5660], 2, doPlot=False, doErr=False, doCov=False)
	y_vph_blu_p, _ = poly_fit_throughput([375, 425, 475], [0.6669, 0.7617, 0.5758], 2, doPlot=False, doErr=False, doCov=False)
	
	#Dichroics - idk which are in which spectrograph, just taking a guess
	#Witness_sample(LLAMAS_dichroic_BG)
	lval_bg, steppt_bg, rval_bg, power_bg, _ = sigmoid_fit_throughput_file(_bgdir+"Witness_sample_BG.txt", doPlot=False, doErr=False)
	
	#Actual transmission_reflection(Witness sample coating lot A)
	lval_sl, steppt_sl, rval_sl, power_sl, _ = sigmoid_fit_throughput_file(_sldir+"Actual_transmission_lot_A.txt", doPlot=False, doErr=False)
	
	#Collimators
	#protoLLAMAS_mirror_MIT-Run-9-2730
	#calculated the average straight from there
	
	#Cameras
	#grab averages from proto files
	y_red_cam = 0.956
	y_blu_cam = 0.954
	#for green, d=0, so I gotta grab the simulated experiment to impute
	gre_priormeans = [prior_mean["gre_l1_t"],prior_mean["gre_l2_t"],prior_mean["gre_l3_t"],prior_mean["gre_l4_t"],prior_mean["gre_l5_t"],prior_mean["gre_l6_t"],prior_mean["gre_l7_t"]]
	y_gre_cam = thru_measurement([], 0, gre_min, gre_max, 0, gre_priormeans, False)

	y_hist_spect1A = [	
			y_gain_red, y_gain_gre, y_gain_blu, 
			y_rn_red, y_rn_gre, y_rn_blu, 
			y_dc_red, y_dc_gre, y_dc_blu, 
			y_qe_red_t[0], y_qe_red_t[1], y_qe_red_t[2], y_qe_red_t[3],	y_qe_red_t[4], 
			y_qe_gre_t[0], y_qe_gre_t[1], y_qe_gre_t[2], y_qe_gre_t[3],	y_qe_gre_t[4], 
			y_qe_blu_t[0], y_qe_blu_t[1], y_qe_blu_t[2], y_qe_blu_t[3],	y_qe_blu_t[4], 
			y_vph_red_p[0], y_vph_red_p[1], y_vph_red_p[2], 
			y_vph_gre_p[0], y_vph_gre_p[1], y_vph_gre_p[2],
			y_vph_blu_p[0], y_vph_blu_p[1], y_vph_blu_p[2], 
			lval_sl, steppt_sl, rval_sl, power_sl, #"y_sl_t0", #"y_sl_t1", #"y_sl_t2", #"y_sl_t3", 
			lval_bg, steppt_bg, rval_bg, power_bg, #"y_bg_t0", #"y_bg_t1", #"y_bg_t2", #"y_bg_t3", 
			0.982, #"y_coll_t",
			y_red_cam, y_gre_cam, y_gre_cam,
			0.077 #"y_frd"
		]
		
	return y_hist_spect1A

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Function to run for this vvopt analysis')
	parser.add_argument('-n', type=int, default=0, help='Number of iterations to give to the function')
	args = parser.parse_args()
	
	###Problem Definition
	d_historical = [
						20,   #t_gain
						30,   #I_gain
						1,	#n_meas_rn
						8,	#d_num
						9600, #d_max
						2,	 #d_pow   #approx, 
						12,  #n_qe #questionable; kind of wasn't measured
						2,  #t_qe
						3, #d_vph_n_pts
						801, #d_dichroic_n_pts
						1501, #d_coll_n_pts #MIT Run 9-2730		
						11, #d_redcam_n_pts #protoLLAMAS camera testing
						0, #d_greencam_n_pts #protoLLAMAS camera testing
						11, #d_bluecam_n_pts #protoLLAMAS camera testing
						10 #d_frd_n_meas #from evaluating_cleaving_through_bigger.xlsx
					]
					
	d_min = [
						0, #t_gain 				0-length exposure for gain exp
						0, #I_gain 				no measurements for gain exp
						0, #n_meas_rn 			no measurements for rn exp
						0, #d_num 				no measurements for dc exp
						0, #d_max 				dc exp filer
						0, #d_pow 				dc exp filler
						0, #n_qe 				no measurements for qe exp
						0, #t_qe 				0-length exposure for qe exp
						0, #d_vph_n_pts 		no measurements for vph exp
						0, #d_dichroic_n_pts 	no measurements for dichroic exp
						0, #d_coll_n_pts		no measurements for collimator exp
						0, #d_redcam_n_pts		no measurements for camera exp
						0, #d_greencam_n_pts	no measurements for camera exp
						0, #d_bluecam_n_pts		no measurements for camera exp
						0  #d_frd_n_meas 		no measurements for FRD exp
					]
	problem = construct_llamas_snr_problem()
	
	req = 3.0
	theta_nominal = problem.theta_nominal
	y_nominal = problem.eta(theta_nominal, d_historical, err=False)

	###Make y_hist
	if args.run == "y_hist":
		y = get_y_hist(problem)
		print(y)

	###Propagate through the joint model to get p(Q|y_hist)
	
	
	###Propagate d* through the exp model to get p(y|theta,d*) and its mean, y*
	
	###Propagate y* through the joint model to get p(Q| y*)
	
	###Compare