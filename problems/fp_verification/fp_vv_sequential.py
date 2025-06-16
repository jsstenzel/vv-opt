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
#focal plane
from problems.fp_verification.fp_problem_sequential import *
#analysis
from obed.obed_multivar import *
from obed.mcmc import *
from obed.obed_gbi import *
from obed.pdf_estimation import *
from uq.uncertainty_propagation import *
from uq.saltelli_gsa import *
from uq.gsa_convergence import *
from uq.gsa_plot import *
#from uq.sensitivity_analysis import *
from inference.bn_modeling import *
from inference.bn_evaluation import *
from opt.utility_max_costcap import *
import pickle
from obed.mcmc import *


################################
#Analysis functions
################################

def vv_nominal(problem, req, theta_nominal, y_nominal):
	print("QoI requirement:", req)
	QoI_nominal = problem.H(theta_nominal, verbose=True)
	print("Given the nominal theta:", theta_nominal)
	#print("Nominal y:", y_nominal)
	print("Nominal QoI:", QoI_nominal)

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Functions to run for this vvopt analysis')
	parser.add_argument('-n', type=int, default=0, help='Number of iterations to give to the function')
	parser.add_argument('-ncomp', type=int, default=0, help='Number of components for GMM')	
	args = parser.parse_args()
	
	###Problem Definition
	problem = fp2
	req = 4.38 #max noise
	d_historical = [
		#20,   #t_gain
		#30,   #I_gain
		1,	  #n_meas_rn
		8,	  #d_num
		9600, #d_max
		2		#d_pow   #approx
	]

	theta_nominal = [2.5, .001]
	QoI_nominal = fp2.H(theta_nominal)
	y_nominal = problem.eta(theta_nominal, d_historical, err=False)

	###Uncertainty Quantification
	if args.run == "nominal":
		vv_nominal(problem, req, theta_nominal, y_nominal)
		
	elif args.run == "UP_jitter_from_BN":
		Qs, _ = bn_load_samples(problem, savefile="BN_sequential_samples", doPrint=True, doDiagnostic=True)
		uncertainty_prop_plot(Qs, xlab="QoI: Avg. Noise [e-]", vline=[req])
		
	elif args.run == "SA_QoI":
		N = args.n
		param_defs = [                             
			["gain", [2.5,0.25**2], "gamma_mv"], #mean, variance
			["rn",   [0.001,.001**2], "gamma_mv"], 
		]
		var_names = [pdef[0] for pdef in param_defs]
		var_bounds = [pdef[1] for pdef in param_defs]
		var_dists = [pdef[2] for pdef in param_defs]

		def model(params):
			theta = params
			return problem.H(theta, verbose=False)
			
		###Generate some new samples and save
		#I want to save samples as often as possible. Therefore, iterate through N two at a time, corresponding to 2(p+2) model evals at a time
		for _ in range(int(N/2)):
			saltelli_eval_sample("SA_FP", 2, var_names, var_dists, var_bounds, model, doPrint=True)
			
	elif args.run == "SA_QoI_evaluate":
		var_names = ["read noise","dark current"]
		
		#S, ST, n_eval = saltelli_indices("SA_QoI", var_names, do_subset=0, doPrint=True)
		total_order_convergence_tests(1200, "SA_FP", var_names, do_subset=0)
		
	elif args.run == "SA_QoI_plot":
		varnames = ["read noise","dark current"]
		S1 = [0.3738, 0.6343]
		ST = [0.3496, 0.6201]
		plot_gsa_full(varnames, S1, ST, S1_conf=[0.0018,0.0013], ST_conf=[0.0040,0.0032], title="", coplot=False, screening=0, xspin=False)
	
	###Train Bayesian network model
	elif args.run == "BN_sample":
		rate = 100
		bn_sampling(problem, savefile="BN_sequential_samples", N=4000000, buffer_rate=rate, doPrint=True)
	
	elif args.run == "BN_train":
		#Train the BN off of the saved data
		ncomp = 45
		q, _ = bn_load_samples(problem, savefile="BN_sequential_samples", doPrint=True, doDiagnostic=True)
		gmm = bn_train_from_file(problem, savefile="BN_sequential_samples", do_subset=args.n, ncomp=ncomp, doPrint=True)
		
		#Save the GMM to a file
		filename = "BN_sequential_model_" + str(len(q)) + "_ncomp" + str(ncomp) + '.pkl'
		#filename = "BN_model.csv"
		bn_save_gmm(gmm, gmm_file=filename)
		
	#Find the highest utility design, subject to a cost cap
	elif args.run == "OPT_costcap":
		print(35482.15)
		minimize_with_penalty(
		problem, 
		costcap=35482.15, 
		gmm_file="BN_sequential_model_4000000_ncomp45.pkl", 
		ylist_file="BN_sequential_samples.csv",
		n_mc=50000, 
		n_tries=1, 
		x0=[14, 21, 4.35234, 2.75273],
		ftol=1e-8, #0.0006407042183632374,
		penalty=10
		)
		
	elif args.run == "prior_update":
		###for the optimal design, d1, assume data that is nominal
		d2 = d_historical #TBD
		#ydata = [yi*1.1 for yi in y_nominal] #test data, just choosing things on the large end of y
		ydata = problem.eta(theta_nominal, d2, err=False)		

		###see what the benefit would be in the batch timeline, i.e. apply the design all at once
		fp2_prior_update(ydata, d_historical, n_mcmc=args.n, loadKDE=True, doDiagnostic=True)

	else:
		print("I don't recognize the command",args.run)
