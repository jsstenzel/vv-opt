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
#Starting a new file just so I can do SA for individual experiment models, and not the whole thing
################################

#sensitivity analysis of HLVA
def vv_SA_rn_red_sample(problem, N=10000):	
	param_defs = [                             
		["gain_red", [0.999,0.2**2], "gamma_mv"],
		["rn_red", [2.32,0.25**2], "gamma_mv"],
		["n_meas", [0,50], "uniform"],
		["sigma_stray", [.001, .01], "uniform"], #for each precision, the mean is the nominal value, and the var is 
		["sigma_dc", [0.3, 0.7], "uniform"]
	]
	var_names = [pdef[0] for pdef in param_defs]
	var_bounds = [pdef[1] for pdef in param_defs]
	var_dists = [pdef[2] for pdef in param_defs]

	def rn_model(params):
		###First, process theta, sampling using the new precisions:
		gain = params[0]
		rn = params[1]
		n_meas = params[2] #this is known precisely
		sigma_stray = params[3]
		sigma_dc = params[4]
		t = 0.1
		
		prior_mean = 2.32 #rn_red
		nx = 2048
		ny = 2048
		
		if n_meas <= 0:
			#This means we don't run the experiment
			#Model the imputation that would occur: 
			# - Ignore the provided theta and assume theta is the mean of the prior
			# - calculate the y that would result in
			return gain * prior_mean

		#Define design variables
		#n_meas number of measurements
		
		sigma_si = gain * math.sqrt(rn**2 + (sigma_dc*t)**2 + (sigma_stray*t)**2)
		
		n = nx * ny * n_meas
		random_var = (2/n) * sigma_si**4
		random_sigma = math.sqrt(random_var)
		random = scipy.stats.norm.rvs(scale = random_sigma)
		y = sigma_si + random
		return y	
		
	###Generate some new samples and save
	#I want to save samples as often as possible. Therefore, iterate through N two at a time, corresponding to 2(p+2) model evals at a time
	for _ in range(int(N/2)):
		saltelli_eval_sample("SA_exp_rn", 2, var_names, var_dists, var_bounds, rn_model, doPrint=True)

#sensitivity analysis
def vv_SA_rn_red_evaluate():
	var_names = ["gain_red","rn_red","n_meas","sigma_stray","sigma_dc"]

	#S, ST, n_eval = saltelli_indices("SA_QoI", var_names, do_subset=0, doPrint=True)
	total_order_convergence_tests(1200, "SA_exp_rn", var_names, do_subset=0)
	
def vv_SA_dc_red_sample(problem, N=10000):	
	param_defs = [                             
		["gain_red", [0.999,0.2**2], "gamma_mv"],
		["rn_red", [2.32,0.25**2], "gamma_mv"],
		["dc_red", [0.00238,.001**2], "gamma_mv"],
		["d_num", [0, 25], "uniform"],      #dc
		["d_max", [1, 12000], "uniform"], #dc
		["d_pow", [0.01,3], "uniform"],      #dc
		["sigma_stray", [.001, .01], "uniform"], #for each precision, the mean is the nominal value, and the var is 
		["sigma_dc", [0.3, 0.7], "uniform"]
	]
	var_names = [pdef[0] for pdef in param_defs]
	var_bounds = [pdef[1] for pdef in param_defs]
	var_dists = [pdef[2] for pdef in param_defs]

	def dc_model(params):
		gain = params[0]
		rn = params[1]
		dc = params[2]
		d_num = math.floor(params[3])
		d_max = math.floor(params[4])
		d_pow = params[5]
		mu_stray = 0
		sigma_stray = params[6]
		sigma_dc = params[7]
		t_0 = 0.1
	
		prior_mean = 0.00238 #dc_red
		nx = 2048
		ny = 2048
		
		if d_num <= 1:
			return prior_mean
		
		#handle design-derived variables
		t_list = []
		for i in range(d_num):
			t = dark_current_time_fn(i, tmin=t_0, dmax=d_max, dpow=d_pow, dnum=d_num)
			t_list.append(t)
			
		#Calculate mean signal at each datapoint
		signal_list = [0]*len(t_list)
		for i,_ in enumerate(signal_list):
			signal_list[i] = gain * (dc + mu_stray) * t_list[i]
		
		#Apply error to each datapoint - you could do it at the end, but why not here.
		for i,_ in enumerate(signal_list):
			ti = t_list[i]
			random_var = gain**2 * (rn**2 + (sigma_dc*ti)**2 + (sigma_stray*ti)**2) / np.sqrt(nx**2 * ny**2)
			random_sigma = math.sqrt(random_var)
			signal_list[i] += scipy.stats.norm.rvs(scale = random_sigma) #apply error to mean
			
		#Calculate slope (assume an intercept but don't worry abot it)
		Sxy, Sx2 = 0,0
		xmean = np.mean(t_list)
		ymean = np.mean(signal_list)
		for i,xi in enumerate(t_list):
			yi = signal_list[i]
			Sxy += (xi - xmean)*(yi - ymean)
			Sx2 += (xi - xmean)**2
		
		m = Sxy / Sx2
		return m
		
	###Generate some new samples and save
	#I want to save samples as often as possible. Therefore, iterate through N two at a time, corresponding to 2(p+2) model evals at a time
	for _ in range(int(N/2)):
		saltelli_eval_sample("SA_exp_dc", 2, var_names, var_dists, var_bounds, dc_model, doPrint=True)

#sensitivity analysis of HLVA
def vv_SA_dc_red_evaluate():
	var_names = ["gain_red","rn_red","dc_red","d_num","d_max","d_pow","sigma_stray","sigma_dc"]

	#S, ST, n_eval = saltelli_indices("SA_QoI", var_names, do_subset=0, doPrint=True)
	total_order_convergence_tests(1200, "SA_exp_dc", var_names, do_subset=0)

def vv_SA_QoI_convergence():
	var_names = ["gain_red","gain_gre","gain_blu","rn_red","rn_gre","rn_blu","dc_red","dc_gre","dc_blu","qe_red_prec","qe_gre_prec","qe_blu_prec","vph_red_prec","vph_gre_prec","vph_blu_prec","sl_prec","bg_prec","coll_prec","red_l1_prec","red_l2_prec","red_l3_prec","red_l4_prec","red_l5_prec","red_l6_prec","red_l7_prec","gre_l1_prec","gre_l2_prec","gre_l3_prec","gre_l4_prec","gre_l5_prec","gre_l6_prec","gre_l7_prec","blu_l1_prec","blu_l2_prec","blu_l3_prec","blu_l4_prec","blu_l5_prec","blu_l6_prec","blu_l7_prec","blu_l8_prec","fiber_frd"]
	
	###Perform the analysis at a few evaluation points
	list_S = []
	list_ST = []
	list_n = [10,50,100,500,1000,5000,10000,50000]
	for n in list_n:
		#S, ST, n_eval = saltelli_indices("SA_QoI", var_names, do_subset=0, doPrint=True)
		#list_S.append(S)
		#list_ST.append(ST)
		total_order_convergence_tests(800, "SA_QoI", var_names, do_subset=n**2)


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
						0, #d_greencam_n_pts		no measurements for camera exp
						0, #d_bluecam_n_pts		no measurements for camera exp
						0  #d_frd_n_meas 		no measurements for FRD exp
					]
	
	###Dodge problem construction
	if args.run == "rn_evaluate":
		vv_SA_rn_red_evaluate()
		sys.exit()
		
	if args.run == "dc_evaluate":
		vv_SA_dc_red_evaluate()
		sys.exit()
		
	problem = construct_llamas_snr_problem()
	
	req = 3.0
	theta_nominal = problem.theta_nominal
	y_nominal = problem.eta(theta_nominal, d_historical, err=False)

	if args.run == "rn_sample":
		vv_SA_rn_red_sample(problem, args.n)
		
	if args.run == "dc_sample":
		vv_SA_dc_red_sample(problem, args.n)