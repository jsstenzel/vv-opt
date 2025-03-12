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

import cProfile
import pstats

def profile_func(func, args, filename, numfunc=30):
	#profile with cProfile
	def function_to_profile():
		return func(*args)
	cProfile.runctx('function_to_profile()', globals(), locals(), filename+'.prof')
	
	#analyze with pstats
	p = pstats.Stats(filename+'.prof')
	p.sort_stats('time').print_stats(numfunc)  # Print top 10 functions by time
	
def bootstrap_profile():
	var_names = ["gain_red","gain_gre","gain_blu","rn_red","rn_gre","rn_blu","dc_red","dc_gre","dc_blu","qe_red_prec","qe_gre_prec","qe_blu_prec","vph_red_prec","vph_gre_prec","vph_blu_prec","sl_prec","bg_prec","coll_prec","red_l1_prec","red_l2_prec","red_l3_prec","red_l4_prec","red_l5_prec","red_l6_prec","red_l7_prec","gre_l1_prec","gre_l2_prec","gre_l3_prec","gre_l4_prec","gre_l5_prec","gre_l6_prec","gre_l7_prec","blu_l1_prec","blu_l2_prec","blu_l3_prec","blu_l4_prec","blu_l5_prec","blu_l6_prec","blu_l7_prec","blu_l8_prec","fiber_frd"]

	S, ST, n_eval = saltelli_indices("SA_QoI", var_names, do_subset=10000, doPrint=False)
	total_order_convergence_tests(100, "SA_QoI", var_names, do_subset=10000)

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Functions to run for this vvopt analysis')
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
						0, #d_greencam_n_pts		no measurements for camera exp
						0, #d_bluecam_n_pts		no measurements for camera exp
						0  #d_frd_n_meas 		no measurements for FRD exp
					]
	problem = construct_llamas_snr_problem()
	
	theta_nominal = problem.theta_nominal
	y_nominal = problem.eta(theta_nominal, d_historical, err=False)
	
	x_dict = deepcopy(problem.x_dict)
	prior_mean_dict = deepcopy(problem.prior_mean_dict)
	
	###Diagnostics
	if args.run == "eta_profile" :
		profile_func(problem.eta, [theta_nominal, d_historical, [], True], 'eta_stats')
		
	elif args.run == "eta_dmin_profile" :
		profile_func(problem.eta, [theta_nominal, d_min, [], True], 'eta_dmin_stats')
		
	elif args.run == "H_profile":
		profile_func(problem.H, [theta_nominal, [], False], 'H_stats')
		
	elif args.run == "ptheta_profile":
		profile_func(problem.prior_rvs, [1], 'H_theta_stats')
		
	elif args.run == "SA_QoI_bootstrap_profile":
		profile_func(bootstrap_profile, [], 'SA_QoI_bootstrap_stats')
		
	elif args.run == "BN_test":
		gmm = bn_load_gmm("BN_model_1639027_ncomp200.pkl")
		presampled_ylist = bn_load_y(problem, "BN_samples_1639027.csv", do_subset=10000, doPrint=False, doDiagnostic=False)
		
		profile_func(U_varH_gbi_joint_presampled, [d_historical, problem, gmm, presampled_ylist, 10000, True], "BN_test_stats")
	
		
	elif args.run == "gain_exp":
		#Check to make sure the function works as expected over a range of inputs
		
		#First, for nominal design and no err, see range of gain						
		t_gain_hist = 20   #t_gain
		I_gain_hist = 30   #I_gain
		
		gains = np.linspace(0, 5, 5000)
		y_gains = []
		y_gains_err = []
		for g in gains:
			yg = gain_exp(g, prior_mean_dict["rn_red"], prior_mean_dict["dc_red"], t_gain_hist, I_gain_hist, x_dict, prior_mean_dict["gain_red"], err=False)
			y_gains.append(yg)
			yge = [gain_exp(g, prior_mean_dict["rn_red"], prior_mean_dict["dc_red"], t_gain_hist, I_gain_hist, x_dict, prior_mean_dict["gain_red"], err=True) for _ in range(100)]
			yge_mean = np.mean(yge)
			y_gains_err.append(yge_mean)
			
		plt.plot(gains, y_gains_err, c='maroon')
		plt.plot(gains, y_gains, c='darkblue')
		plt.legend(["err=True","err=False"])
		plt.xlabel("theta_gain")
		plt.ylabel("y_gain")
		plt.title("Output of gain_exp for red camera, historical design")
		plt.show()
		
		plt.plot(gains, [abs(ye - y) for ye,y in zip(y_gains_err,y_gains)], c='red')
		plt.xlabel("theta_gain")
		plt.ylabel("Mean absolute error")
		plt.title("Output of gain_exp for red camera, historical design")
		plt.show()
		
		#Now, show 2d heatmap of how the design variables impact nominal/historical
		ts = np.linspace(.1, 600, 100)
		Is = np.linspace(1, 100, 99)
		
		y_gains = np.empty((100,99))
		for it, t in enumerate(ts):
			for iI, I in enumerate(Is):
				#yg = gain_exp(g, prior_mean_dict["rn_red"], prior_mean_dict["dc_red"], t_gain_hist, I_gain_hist, x_dict, prior_mean_dict["gain_red"], err=False)
				yge = [gain_exp(g, prior_mean_dict["rn_red"], prior_mean_dict["dc_red"], t_gain_hist, I_gain_hist, x_dict, prior_mean_dict["gain_red"], err=True) for _ in range(1000)]
				y_gains[it,iI] = np.var(yge, ddof=1)
		sns.heatmap(y_gains, cmap="inferno")
		plt.title("Variance of gain_exp function at 100 samples")
		plt.ylabel('t')
		plt.xlabel('I')
		plt.show()
		
	else:
		print("Profile method not recognized.")