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

def vv_nominal(problem, req, theta_nominal, y_nominal):
	print(problem)
	print("QoI requirement:", req)
	QoI_nominal = problem.H(theta_nominal, verbose=True)
	print("Given the nominal theta:", theta_nominal)
	print("Nominal y:", y_nominal)
	print("Nominal QoI:", QoI_nominal)

###uncertainty analysis
def vv_UA_theta(problem, n=10**4):
	#uncertainty propagation of HLVA
	uq_thetas = problem.prior_rvs(n)
	
	covmatrix_heatmap(uq_thetas, names=problem.theta_names, rescale=True)
	
def vv_UP_QoI(problem, req, n=10**4):
	#uncertainty propagation of HLVA
	uq_thetas = problem.prior_rvs(n)
	Qs = []
	for i,theta in enumerate(uq_thetas):
		print(i, flush=True)
		try:
			Q = problem.H(theta) 
			print(Q, flush=True)
			Qs.append(Q)
		except:
			print("System model eval fail?", flush=True)
	print(Qs)
	#uncertainty_prop_plot([theta[0] for theta in uq_thetas], xlab="How to plot this...")
	uncertainty_prop_plot(Qs, xlab="QoI: SNR", vline=[req])

	#prob of meeting req along priors:
	count_meetreq = 0
	for Q in Qs:
		if Q >= req:
			count_meetreq += 1
	prob_meetreq = count_meetreq / len(Qs)
	print("Probability of meeting requirement given priors:", prob_meetreq)
	
def vv_UP_QoI_samples(req, base_name="SA_QoI", doPrint=True, do_subset=0):
	var_names = ["gain_red","gain_gre","gain_blu","rn_red","rn_gre","rn_blu","dc_red","dc_gre","dc_blu","qe_red_prec","qe_gre_prec","qe_blu_prec","vph_red_prec","vph_gre_prec","vph_blu_prec","sl_prec","bg_prec","coll_prec","red_l1_prec","red_l2_prec","red_l3_prec","red_l4_prec","red_l5_prec","red_l6_prec","red_l7_prec","gre_l1_prec","gre_l2_prec","gre_l3_prec","gre_l4_prec","gre_l5_prec","gre_l6_prec","gre_l7_prec","blu_l1_prec","blu_l2_prec","blu_l3_prec","blu_l4_prec","blu_l5_prec","blu_l6_prec","blu_l7_prec","blu_l8_prec","fiber_frd"]

	###Make sure the files exist
	if not os.path.isfile(base_name+'_A.csv'):
		print("File",base_name+'_A.csv',"is missing")
		sys.exit()
	if not os.path.isfile(base_name+'_B.csv'):
		print("File",base_name+'_B.csv',"is missing")
		sys.exit()
	
	###Safely read out all of the samples into matrices
	Ay = []
	By = []
	
	if doPrint:
		print("Reading the data files...",flush=True)
	
	if do_subset == 0:
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for l,row in enumerate(csvreader):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_A.csv')
				else:
					Ay.append([float(elem) for elem in row])
		with open(base_name+'_B.csv') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for l,row in enumerate(csvreader):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_B.csv')
				else:
					By.append([float(elem) for elem in row])
	else:
		lim = int(do_subset/2)
		###Optionally, we can analyze less than the full set of provided samples
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for l,row in enumerate(islice(csvreader, lim)):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_A.csv')
				else:
					Ay.append([float(elem) for elem in row])
		
		with open(base_name+'_B.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for l,row in enumerate(islice(csvreader, lim)):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_B.csv')
				else:
					By.append([float(elem) for elem in row])
	
	#Pull out all of the y samples
	Ay.extend(By)
	Q_samples = [Ay_row[-1] for Ay_row in Ay] #only last element
	
	#Do the UQ
	uncertainty_prop_plot(Q_samples, xlab="QoI: SNR", vline=[req])

	#prob of meeting req along priors:
	count_meetreq = 0
	for Q in Q_samples:
		if Q >= req:
			count_meetreq += 1
	prob_meetreq = count_meetreq / len(Q_samples)
	print("Probability of meeting requirement given priors:", prob_meetreq)

#sensitivity analysis of HLVA
def vv_SA_QoI_sample(problem, N=10000):
	#Set up the problem. The parameters are mostly theta,
	#Except replace every Gaussian Process with a prior with a hyperparameter on the precision, p=1/sigma^2
	#And put a gamma distribution as the prior for that, according to the principle of maximum entropy
	#Note that all of these variances are in u-space
	sloc = './priors/'
	gp_variance_list = [
		list(load_gp_prior_from_file(sloc+'prior_gp_qe_red'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gp_qe_gre'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gp_qe_blu'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gp_vph_red'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gp_vph_gre'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gp_vph_blu'))[0],
		1.0, #sl,
		1.0, #bg
		1.0, #coll
		list(load_gp_prior_from_file(sloc+'prior_red1'))[0],
		list(load_gp_prior_from_file(sloc+'prior_red2'))[0],
		list(load_gp_prior_from_file(sloc+'prior_red3'))[0],
		list(load_gp_prior_from_file(sloc+'prior_red4'))[0],
		list(load_gp_prior_from_file(sloc+'prior_red5'))[0],
		list(load_gp_prior_from_file(sloc+'prior_red6'))[0],
		list(load_gp_prior_from_file(sloc+'prior_red7'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gre1'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gre2'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gre3'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gre4'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gre5'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gre6'))[0],
		list(load_gp_prior_from_file(sloc+'prior_gre7'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu1'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu2'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu3'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu4'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu5'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu6'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu7'))[0],
		list(load_gp_prior_from_file(sloc+'prior_blu8'))[0]
	]
	gp_prec_list = [1/v for v in gp_variance_list]
	hyper_mean, hyper_std = uncertainty_prop(gp_prec_list, False, False)
	#This is the best info I have about the priors, so I should use it in the development of a prior for the hyperparameters
	
	param_defs = [                             
		["gain_red", [0.999,0.2**2], "gamma_mv"],
		["gain_gre", [1.008,0.2**2], "gamma_mv"],
		["gain_blu", [1.008,0.2**2], "gamma_mv"],
		["rn_red", [2.32,0.25**2], "gamma_mv"],
		["rn_gre", [2.35,0.25**2], "gamma_mv"],
		["rn_blu", [2.35,0.25**2], "gamma_mv"],
		["dc_red", [0.00238,.001**2], "gamma_mv"],
		["dc_gre", [0.00267,.001**2], "gamma_mv"],
		["dc_blu", [0.00267,.001**2], "gamma_mv"],
		["qe_red_prec", [hyper_mean, hyper_std**2], "gamma_mv"], #for each precision, the mean is the nominal value, and the var is 
		["qe_gre_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["qe_blu_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["vph_red_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["vph_gre_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["vph_blu_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["sl_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["bg_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["coll_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["red_l1_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["red_l2_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["red_l3_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["red_l4_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["red_l5_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["red_l6_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["red_l7_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["gre_l1_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["gre_l2_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["gre_l3_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["gre_l4_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["gre_l5_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["gre_l6_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["gre_l7_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l1_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l2_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l3_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l4_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l5_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l6_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l7_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["blu_l8_prec", [hyper_mean, hyper_std**2], "gamma_mv"],
		["fiber_frd", [3.434339623321249, 133.9392453095287], "beta"]
	]
	var_names = [pdef[0] for pdef in param_defs]
	var_bounds = [pdef[1] for pdef in param_defs]
	var_dists = [pdef[2] for pdef in param_defs]

	def model(params):
		###First, process theta, sampling using the new precisions:
		theta = []
		for i,name in enumerate(var_names):
			if name.endswith('_prec'):
				variance = params[i] #the hyperparameter
				ls = problem.priors[i][1][1]
				prior_pts = problem.priors[i][1][2]
				mean_fn = problem.priors[i][1][3]

				gp_prior = GaussianProcessDist1D(variance, ls, prior_pts, mean_fn)
				sample = gp_prior.sample()
				theta.append(sample)
			else:
				theta.append(params[i])
		
		###Finally, return QoI
		return problem.H(theta, verbose=False)
		
	###Generate some new samples and save
	#I want to save samples as often as possible. Therefore, iterate through N two at a time, corresponding to 2(p+2) model evals at a time
	for _ in range(int(N/2)):
		saltelli_eval_sample("SA_QoI", 2, var_names, var_dists, var_bounds, model, doPrint=True)

#sensitivity analysis of HLVA
def vv_SA_QoI_evaluate(problem):
	var_names = ["gain_red","gain_gre","gain_blu","rn_red","rn_gre","rn_blu","dc_red","dc_gre","dc_blu","qe_red_prec","qe_gre_prec","qe_blu_prec","vph_red_prec","vph_gre_prec","vph_blu_prec","sl_prec","bg_prec","coll_prec","red_l1_prec","red_l2_prec","red_l3_prec","red_l4_prec","red_l5_prec","red_l6_prec","red_l7_prec","gre_l1_prec","gre_l2_prec","gre_l3_prec","gre_l4_prec","gre_l5_prec","gre_l6_prec","gre_l7_prec","blu_l1_prec","blu_l2_prec","blu_l3_prec","blu_l4_prec","blu_l5_prec","blu_l6_prec","blu_l7_prec","blu_l8_prec","fiber_frd"]

	total_order_convergence_tests(1200, "SA_QoI", var_names, do_subset=0)
	sys.exit()

	###Perform the analysis at a few evaluation points
	list_S = []
	list_ST = []
	list_n = [10,100,1000,10000]
	for n in list_n:
		S, ST, n_eval = saltelli_indices("SA_QoI", var_names, do_subset=0, doPrint=True)
		list_S.append(S)
		list_ST.append(ST)

#Uncertainty analysis of the experiment models
def vv_UP_exp(problem, dd, theta_nominal, n=10**4, savefig=False):
	print("Likelihood distribution for nominal historical case:",flush=True)
	#tt = problem.prior_rvs(1); print(tt)
	tt = theta_nominal
	ysample_nominal = [problem.eta(tt, dd) for _ in range(n)]
	uncertainty_prop_plots(ysample_nominal, xlabs=problem.y_names, saveFig='UP_exp_nominal' if savefig else '')
	
	#Also plot the y's generated by the joint distribution p(y|theta,d)p(theta)
	print("Experiments simulated from the joint distribution:",flush=True)
	uq_thetas = problem.prior_rvs(n)
	uq_ys = [problem.eta(theta, dd) for theta in uq_thetas]
	uncertainty_prop_plots(uq_ys, c='orchid', xlabs=problem.y_names, saveFig='UP_joint_y0' if savefig else '')

#These cause problems on eofe8
#sensitivity analysis of the experiment models
def vv_SA_exp(problem, dd, p=8):
	#we want to see sensitivity of each yi to their inputs, theta and x
	#challenge: sobol_saltelli expects a function that takes a single list of parameters
	#right now, im doing this in a way that just deals with the thetas
	#i think i can handle including x's intelligently as well some day, using this filter:
	SA_filter = range(len(problem.theta_names)) #[]
		
	for i,yi in enumerate(problem.y_names):
		def exp_fn_i(param):
			y = problem.eta(param, dd) #want to figure out a smarter way to do analysis of x in the future
			yi = y[i]
			return yi

		expvar_names = problem.theta_names #+ problem.x_names
		expvar_dists = [prior[0] for prior in problem.priors] #+ [prior[0] for prior in problem.x_dists]
		expvar_bounds=[prior[1] for prior in problem.priors] #+ [prior[1] for prior in problem.x_dists]
			
		#so ugly!!!
		Si_1 = sobol_saltelli(exp_fn_i, 
							2**p, #SALib wants powers of 2 for convergence
							var_names=[x for i,x in enumerate(expvar_names) if i in SA_filter], 
							var_dists=[x for i,x in enumerate(expvar_dists) if i in SA_filter], 
							var_bounds=[x for i,x in enumerate(expvar_bounds) if i in SA_filter], 
							conf = 0.99, doSijCalc=False, doPlot=True, doPrint=True)	
							
def vv_SA_joint(problem, p=5):
	#it'll be straightforward to see the dependence of QoI on theta
	Si = sobol_saltelli(problem.H, 
						2**p, #SALib wants powers of 2 for convergence
						var_names=problem.theta_names, 
						var_dists=[prior[0] for prior in problem.priors], 
						var_bounds=[prior[1] for prior in problem.priors], #need to do something odd here....?
						conf = 0.95, doSijCalc=False, doPlot=True, doPrint=True)

def vv_gbi_test(problem, d, N, y=[], ncomp=0):		
	print("Training...")
	theta_train = problem.prior_rvs(N)
	qoi_train = [problem.H(theta) for theta in theta_train]
	y_train = [problem.eta(theta, d) for theta in theta_train]
	
	gmm = gbi_train_model(theta_train, qoi_train, y_train, verbose=2, ncomp=ncomp)
	
	if y==[]:
		print("Determining random measurement Yd...")
		truth_theta = problem.prior_rvs(1)
		Yd = problem.eta(truth_theta, d)
	else:
		Yd = y
	
	print("Conditioning...")
	a,b,c = gbi_condition_model(gmm, Yd, verbose=2)
	
	if y==[]:
		plot_predictive_posterior(a, b, c, 0, 500, drawplot=False, plotmean=True)
		plt.axvline(problem.H(truth_theta), c='blue')
		plt.show()
	else:
		plot_predictive_posterior(a, b, c, 0, 500, drawplot=True)

def vv_obed_gbi(problem, d):
	U, U_list = U_varH_gbi(d, problem, n_mc=10**4, n_gmm=10**4, doPrint=True)
	print(U)
	#print(U_list)
	uncertainty_prop_plot(U_list, c='royalblue', xlab="specific U")#, saveFig='OBEDresult')
	return U
	
def uncertainty_mc(problem, dd, n_mc=10**2, n_gmm=10**2, n_test=10**2):
	util_samples = []
	for ii in range(n_test):
		print(ii, flush=True)
		util, _ = U_varH_gbi(dd, problem, n_mc=n_mc, n_gmm=n_gmm, ncomp=10, doPrint=False)
		util_samples.append(util)
		
	uncertainty_prop_plot(util_samples, c='purple', xlab="utility for d=d_hist")
	print(statistics.variance(util_samples))
	return util_samples

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
	problem = construct_llamas_snr_problem()
	
	req = 3.0
	theta_nominal = problem.theta_nominal
	y_nominal = problem.eta(theta_nominal, d_historical, err=False)

	###Uncertainty Quantification
	if args.run == "nominal":
		vv_nominal(problem, req, theta_nominal, y_nominal)
		
	if args.run == "system_model":
		if args.n <= 1:
			theta = problem.prior_rvs(1)
			problem.print_theta(theta)
			QoI = problem.H(theta, verbose=True)
			print("QoI:",QoI,flush=True)
		else:
			thetas = problem.prior_rvs(args.n)
			for t in thetas:
				QoI = problem.H(t, verbose=True)
				print("QoI:",QoI,flush=True)
		
	if args.run == "cheap_design":
		print("Given the nominal theta:", theta_nominal)
		print("and the cheapest design:", d_min)
		y_cheap = problem.eta(theta_nominal, d_min, err=True)
		
		print("Cheapest y:", y_cheap)
		print("Cost of design:", problem.G(d_min))
	
	if args.run == "UA_theta":
		vv_UA_theta(problem, n=args.n)
	
	if args.run == "UP_QoI":
		vv_UP_QoI(problem, req, n=args.n)
		
	if args.run == "UP_QoI_samples":
		vv_UP_QoI_samples(req, base_name="SA_QoI", doPrint=True)
		
	if args.run == "UP_exp":
		vv_UP_exp(problem, d_historical, theta_nominal, n=args.n)
	
	if args.run == "SA_QoI_sample":
		vv_SA_QoI_sample(problem, N=args.n)

	if args.run == "SA_QoI_evaluate":
		vv_SA_QoI_evaluate(problem)
	
	#Still needs massaging...
	#if args.run == "SA_exp":
	#	vv_SA_exp(problem, d_historical)

	###Optimal Bayesian Experimental Design
	if args.run == "BN_sample":
		bn_sampling(problem, savefile="BN_samples", N=args.n, doPrint=True)
	
	if args.run == "BN_train":
		#Train the BN off of the saved data
		q, _ = bn_load_samples(problem, savefile="BN_samples", doPrint=True, doDiagnostic=True)
		gmm = bn_train_from_file(problem, savefile="BN_samples", doPrint=True)
		
		#Run the validation test
		bn_measure_model_mse(problem, gmm, N=args.n, doPrint=True)
		
		#Save the GMM to a file
		#filename = "BN_" + str(len(q)) + '.csv'
		filename = "BN_model.csv"
		bn_save_gmm(gmm, gmm_file=filename)
	
	if args.run == "OBED_test":
		#Load the GMM from file
		gmm = bn_load_gmm("BN_model.csv")
	
		#Calculate U for several different designs
		U_varH_gbi_joint(d_historical, problem, gmm, n_mc=args.n, ncomp=0, doPrint=True)
		U_varH_gbi_joint(d_min, problem, gmm, n_mc=args.n, ncomp=0, doPrint=True)
	
	"""
	if args.run == "gbi_test":
		vv_gbi_test(problem, d_historical, 10**1, y_nominal, ncomp=0)
	if args.run == "gbi_test_rand":
		vv_gbi_test(problem, d_historical, 10**1, ncomp=10)
	
	if args.run == "obed_gbi":
		U_hist = vv_obed_gbi(problem, d_historical)
	
	if args.run == "uncertainty_mc":
		util_samples = uncertainty_mc(problem, d_historical, n_mc=10**2, n_gmm=10**2, n_test=3)
		print(util_samples)
		
	if args.run == "vv_opt_parallel":
		costs, utilities, designs = ngsa2_problem_parallel(8, problem, hours=0, minutes=0, popSize=10, nMonteCarlo=10**3, nGMM=10**3)
		plot_ngsa2(costs, utilities, showPlot=True, savePlot=False, logPlotXY=[False,False])
	"""
