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
from problems.problem_definition import *
from problems.llamas_snr_full.snr_cost_model import *
#analysis
from obed.obed_gbi import *
#from obed.pdf_estimation import *
from inference.bn_modeling import *
from inference.bn_evaluation import *
from uq.uncertainty_propagation import *
from opt.nsga import *

#This version of the optimization code does not include or pass around the problem object
#this is intended to cut down on pickling time & memory incurred by multiprocessing

def construct_slim_llamas_snr_problem():
	empty_prior = ["uniform", [0,1]]
	theta_defs = [                             #mean, variance
						["gain_red", empty_prior, "continuous"],
						["gain_gre", empty_prior, "continuous"],
						["gain_blu", empty_prior, "continuous"],
						["rn_red", empty_prior, "continuous"],
						["rn_gre", empty_prior, "continuous"],
						["rn_blu", empty_prior, "continuous"],
						["dc_red", empty_prior, "continuous"],
						["dc_gre", empty_prior, "continuous"],
						["dc_blu", empty_prior, "continuous"],
						["qe_red_t", empty_prior, "continuous"],
						["qe_gre_t", empty_prior, "continuous"],
						["qe_blu_t", empty_prior, "continuous"],
						["vph_red_t", empty_prior, "continuous"],
						["vph_gre_t", empty_prior, "continuous"],
						["vph_blu_t", empty_prior, "continuous"],
						["sl_t", empty_prior, "continuous"],
						["bg_t", empty_prior, "continuous"],
						["coll_t", empty_prior, "continuous"],
						["red_l1_t", empty_prior, "continuous"],
						["red_l2_t", empty_prior, "continuous"],
						["red_l3_t", empty_prior, "continuous"],
						["red_l4_t", empty_prior, "continuous"],
						["red_l5_t", empty_prior, "continuous"],
						["red_l6_t", empty_prior, "continuous"],
						["red_l7_t", empty_prior, "continuous"],
						["gre_l1_t", empty_prior, "continuous"],
						["gre_l2_t", empty_prior, "continuous"],
						["gre_l3_t", empty_prior, "continuous"],
						["gre_l4_t", empty_prior, "continuous"],
						["gre_l5_t", empty_prior, "continuous"],
						["gre_l6_t", empty_prior, "continuous"],
						["gre_l7_t", empty_prior, "continuous"],
						["blu_l1_t", empty_prior, "continuous"],
						["blu_l2_t", empty_prior, "continuous"],
						["blu_l3_t", empty_prior, "continuous"],
						["blu_l4_t", empty_prior, "continuous"],
						["blu_l5_t", empty_prior, "continuous"],
						["blu_l6_t", empty_prior, "continuous"],
						["blu_l7_t", empty_prior, "continuous"],
						["blu_l8_t", empty_prior, "continuous"],
						["fiber_frd", empty_prior, "continuous"]
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
	
	x_defs = [
				#rn exp
				["t_rn", [], "continuous", .1], #100 ms exposure
				#dc exp
				["t_0", [], "continuous", 0.1], #100ms baseline exposure assumed
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
	
	eta = None
	H = None
	Gamma = snr_cost
	llamas_snr = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs_discrete, x_defs)
	print("llamas_snr_full SLIM problem constructed.",flush=True)
	return llamas_snr
	
	
def vv_OPT(problem, gmm_file, ysamples_file, design_pts, epsilon, util_err, do_hrs, do_min, threads, popSize, nMC, displayFreq=10):
	#Load the GMM and presampled y from file
	print("Loading GMM and presamples...",flush=True)
	gmm20 = bn_load_gmm(gmm_file)
	presampled_ylist = bn_load_y(problem, ysamples_file, doPrint=False, doDiagnostic=False)
	
	costs, utilities, designs = nsga2_obed_bn(	
					n_threads=threads, 
					prob=problem, 
					hours=do_hrs, 
					minutes=do_min, 
					popSize=popSize, 
					nSkip=2, 
					tolDelta=epsilon,#utility_conf95/1.96, 
					nPeriod=5, 
					nMonteCarlo=nMC, 
					GMM=gmm20, 
					Ylist=presampled_ylist,
					displayFreq=displayFreq
				)
	plot_nsga2(costs, utilities, design_pts, util_err=util_err, showPlot=True, savePlot=False, logPlotXY=[False,False])

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Function to run for this vvopt analysis')
	parser.add_argument('-n', type=int, default=0, help='Number of iterations to give to the function')
	parser.add_argument('--filename', metavar='string', default="SA_QoI", help='Base name to five to SA_QoI_sample')
	args = parser.parse_args()
	
	###Problem Definition
	_wave_min = 350.0
	_wave_max = 975.0
	_wave_redgreen = 690.0
	_wave_greenblue = 480.0
	_bandpass = _wave_max - _wave_min

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
						13, #d_redcam_n_pts #protoLLAMAS camera testing
						0, #d_greencam_n_pts #protoLLAMAS camera testing
						13, #d_bluecam_n_pts #protoLLAMAS camera testing
						10 #d_frd_n_meas #from evaluating_cleaving_through_bigger.xlsx
					]
	d_0 = [
						0.1, #t_gain 				0-length exposure for gain exp
						0, #I_gain 				no measurements for gain exp
						0, #n_meas_rn 			no measurements for rn exp
						0, #d_num 				no measurements for dc exp
						1, #d_max 				dc exp filer
						0.01, #d_pow 				dc exp filler
						0, #n_qe 				no measurements for qe exp
						0.1, #t_qe 				0-length exposure for qe exp
						0, #d_vph_n_pts 		no measurements for vph exp
						0, #d_dichroic_n_pts 	no measurements for dichroic exp
						0, #d_coll_n_pts		no measurements for collimator exp
						0, #d_redcam_n_pts		no measurements for camera exp
						0, #d_greencam_n_pts		no measurements for camera exp
						0, #d_bluecam_n_pts		no measurements for camera exp
						0  #d_frd_n_meas 		no measurements for FRD exp
					]
	d_min = [
						0.1, #t_gain 				0-length exposure for gain exp
						1, #I_gain 				no measurements for gain exp
						1, #n_meas_rn 			no measurements for rn exp
						1, #d_num 				no measurements for dc exp
						1, #d_max 				dc exp filer
						0.01, #d_pow 				more short experiments
						1, #n_qe 				no measurements for qe exp
						0.1, #t_qe 				0-length exposure for qe exp
						1, #d_vph_n_pts 		no measurements for vph exp
						1, #d_dichroic_n_pts 	no measurements for dichroic exp
						1, #d_coll_n_pts		no measurements for collimator exp
						1, #d_redcam_n_pts		no measurements for camera exp
						1, #d_greencam_n_pts		no measurements for camera exp
						1, #d_bluecam_n_pts		no measurements for camera exp
						1  #d_frd_n_meas 		no measurements for FRD exp
					]

	d_max = [
						600, #t_gain
						100, #I_gain
						50, #n_meas_rn
						25, #d_num
						12000, #d_max
						3, #d_pow
						100, #n_qe
						600, #t_qe
						_bandpass*10, #d_vph_n_pts
						_bandpass*10, #d_dichroic_n_pts
						_bandpass*10, #d_coll_n_pts
						(_wave_max-_wave_redgreen)*10, #d_redcam_n_pts
						(_wave_redgreen-_wave_greenblue)*10, #d_greencam_n_pts
						(_wave_greenblue-_wave_min)*10, #d_bluecam_n_pts
						2400  #d_frd_n_meas
					]
	d_med = [dd/2 for dd in d_max]
	d_med[5] = 1 #d_pow
	
	problem = construct_slim_llamas_snr_problem()
	#instead of constructing the problem, im just constructing what nsga needs up top
	
	design_pts = [
		[problem.G(d_historical), 0.004240541527302059, "d_hist", 5.384503718405341e-05],
		[problem.G(d_0), 0.0047856764435419575, "d_0", 5.384503718405341e-05],
		[problem.G(d_max), 0.0022957744137691916, "d_max", 5.384503718405341e-05],
		[problem.G(d_med), 0.004294863943612242, "d_med", 5.384503718405341e-05],
		#[problem.G(d_min), 0.004772754162991483, "d_min", 5.384503718405341e-05]
	]
	
	if args.run == "OPT_test":
		vv_OPT(
			problem,
			gmm_file="ncomp_testing/BN_model_1639027_ncomp20.pkl", 
			ysamples_file="BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=0.01, #probably the ncomp20 model is better than this
			util_err=0.001,
			do_hrs = 0,
			do_min = 1,
			threads = 1 if args.n==0 else args.n,
			popSize=10,# if args.n==0 else args.n,
			nMC=10,
			displayFreq=10
		)

	elif args.run == "OPT_nmc_p4":
		conf95 = 0.0006962967583258008
		conf_frac = conf95 = conf95 / (0.0047856764435419575 - 0.0022957744137691916)
		vv_OPT(
			problem,
			gmm_file="ncomp_testing/BN_model_1639027_ncomp200.pkl", 
			ysamples_file="BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=conf_frac,
			util_err=conf95,
			do_hrs = 11 if args.filename == "timed" else 0,
			do_min = 30 if args.filename == "timed" else 0,
			threads = 5,#1 if args.n == 0 else args.n,
			popSize=50,#30 if args.n==0 else args.n,
			nMC=10**4,
			displayFreq=10
		)

	elif args.run == "OPT_nmc_p5":
		conf95 = 0.00019389166166701476
		conf_frac = conf95 = conf95 / (0.0047856764435419575 - 0.0022957744137691916)
		vv_OPT(
			problem,
			gmm_file="BN_model_1639027_ncomp200.pkl", 
			ysamples_file="ncomp_testing/BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=conf_frac, #rough analogue for nMC=10^5
			util_err=conf95,
			do_hrs = 11 if args.filename == "timed" else 0,
			do_min = 30 if args.filename == "timed" else 0,
			threads = 8 if args.n == 0 else args.n,
			popSize=40,
			nMC=10**5
		)
		
	elif args.run == "OPT_nmc_p6":
		#prepare initializaiton samples based on previous run + prepared designs
	    
		conf95 = 5.384503718405341e-05
		conf_frac = conf95 / (0.0047856764435419575 - 0.0022957744137691916)
		vv_OPT(
			problem,
			gmm_file="BN_model_1639027_ncomp200.pkl", 
			ysamples_file="ncomp_testing/BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=conf_frac, #the result for half the 95% confidence interval width for nMC=10^6
			util_err=conf95,
			do_hrs = 11 if args.filename == "timed" else 0,
			do_min = 30 if args.filename == "timed" else 0,
			threads = 8 if args.n == 0 else args.n,
			popSize=60,
			nMC=10**6,
			samples=[]
		)
	
	else:
		print("I dont recognize the command",args.run)
