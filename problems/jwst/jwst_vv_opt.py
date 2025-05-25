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
#jwst
from problems.jwst.jwst_problem import *
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

def vv_system(problem, theta_nominal):
	print(problem)
	#print("QoI requirement:", req)
	QoI_nominal = problem.H(theta_nominal, verbose=True)
	print("Given the nominal theta:", theta_nominal)
	#print("Nominal y:", y_nominal)
	print("Nominal QoI:", QoI_nominal)
	
def vv_experiment(problem, theta_nominal, d):
	print(problem)
	print("Given the nominal theta:", theta_nominal)
	print("And this design:", theta_nominal)
	ynom = problem.eta(theta_nominal, d, err=False)
	print("Nominal y:", ynom)
	y = problem.eta(theta_nominal, d, err=True)
	print("A sampled y:", y)

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
		print(i)
		try:
			Q = problem.H(theta) 
			print(Q, flush=True)
			Qs.append(Q)
		except:
			print("System model eval fail?", flush=True)
	#print(Qs)
	#uncertainty_prop_plot([theta[0] for theta in uq_thetas], xlab="How to plot this...")
	uncertainty_prop(Qs, xlab="QoI: SNR", vline=[req])

	#prob of meeting req along priors:
	count_meetreq = 0
	for Q in Qs:
		if Q >= req:
			count_meetreq += 1
	prob_meetreq = count_meetreq / len(Qs)
	print("Probability of meeting requirement given priors:", prob_meetreq)


def vv_UP_QoI_samples(base_name="SA_jitter", doPrint=True, do_subset=0):
	var_names = ["Us","Ud","Qc","I_SMhubt","I_SMhuba","K_yPM","I_xRWA","I_yRWA","I_RWt","I_RWa","I_ISOa","I_ISOt","K_yISO","K_xISO","I_bus","I_propt","I_propa","I_i1","I_i2","I_i3","A_sptop","D_sp","t_sp","I_ss","K_rad1","K_rad2","K_rISO","K_act1","K_act2","I_iso","K_zpet","K_pm1","K_pm3","K_pm4","K_pm5","K_pm6","K_act_pm2","K_act_pm3","K_act_pm4","K_act_pm5","K_act_pm6","K_xpet","c_RWA","c_RWAI","c_SM_act","c_PM","c_PM_act","c_petal","zeta_sunshield","zeta_isolator","zeta_solarpanel","Ru","fc","Tst","Srg","Sst","Tgs","lambda_","Ro","QE","Mgs","fca","Kc","Kcf"]

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
	uncertainty_prop(Q_samples, xlab="QoI: Jitter")#, vline=[req])
	
	for i,name in enumerate(var_names):
		t_samples = [Ay_row[i] for Ay_row in Ay]
		uncertainty_prop(t_samples, xlab="UP "+str(name))

	#prob of meeting req along priors:
	#count_meetreq = 0
	#for Q in Q_samples:
	#	if Q >= req:
	#		count_meetreq += 1
	#prob_meetreq = count_meetreq / len(Q_samples)
	#print("Probability of meeting requirement given priors:", prob_meetreq)



#sensitivity analysis of HLVA
def vv_SA_jitter_sample(problem, filename, N=10000):
	#Set up the problem. The parameters are theta and a subset of x

	param_defs = [      
		["Us", unif_margin(1.8), "uniform"],
		["Ud", unif_margin(60.0), "uniform"],
		["Qc", unif_margin(0.005), "uniform"],
		["I_SMhubt", unif_margin(0.25200E+00), "uniform"],
		["I_SMhuba", unif_margin(0.45900E+00), "uniform"],
		["K_yPM", unif_margin(0.77400E+06), "uniform"],
		["I_xRWA", unif_margin(0.40187E+00), "uniform"],
		["I_yRWA", unif_margin(0.22445E+00), "uniform"],
		["I_RWt", unif_margin(0.83595E-01), "uniform"],
		["I_RWa", unif_margin(0.14339E+00), "uniform"],
		["I_ISOa", unif_margin(0.11720E-03), "uniform"],
		["I_ISOt", unif_margin(0.39300E-01), "uniform"],
		["K_yISO", unif_margin(0.14600E+04), "uniform"],
		["K_xISO", unif_margin(0.14000E+12), "uniform"],
		["I_bus", unif_margin(0.85080E+02), "uniform"],
		["I_propt", unif_margin(0.51100E+01), "uniform"],
		["I_propa", unif_margin(0.74000E+00), "uniform"],
		["I_i1", unif_margin(0.49200E+01), "uniform"],
		["I_i2", unif_margin(0.75420E+01), "uniform"],
		["I_i3", unif_margin(0.41280E+01), "uniform"],
		["A_sptop", unif_margin(0.14040E-02), "uniform"],
		["D_sp", unif_margin(0.060), "uniform"],
		["t_sp", unif_margin(0.003), "uniform"],
		["I_ss", unif_margin(0.78350E-08), "uniform"],
		["K_rad1", unif_margin(0.50000E+6), "uniform"],
		["K_rad2", unif_margin(0.30000E+6), "uniform"],
		["K_rISO", unif_margin(3000), "uniform"],
		["K_act1", unif_margin(0.20000E+11), "uniform"],
		["K_act2", unif_margin(0.14000E+12), "uniform"],
		["I_iso", unif_margin(1.00000E-5), "uniform"],
		["K_zpet", unif_margin(0.9000E+08), "uniform"],
		["K_pm1", unif_margin(0.10000E+07), "uniform"],
		["K_pm3", unif_margin(0.58400E+06), "uniform"],
		["K_pm4", unif_margin(0.59820E+02), "uniform"],
		["K_pm5", unif_margin(0.49000E+02), "uniform"],
		["K_pm6", unif_margin(0.33250E+02), "uniform"],
		["K_act_pm2", unif_margin(0.29100E+07), "uniform"],
		["K_act_pm3", unif_margin(0.10000E+07), "uniform"],
		["K_act_pm4", unif_margin(0.33250E+02), "uniform"],
		["K_act_pm5", unif_margin(0.49000E+02), "uniform"],
		["K_act_pm6", unif_margin(0.12012E+03), "uniform"],
		["K_xpet", unif_margin(1e16), "uniform"],
		["c_RWA", unif_margin(4.23245e-7*0.01*math.sqrt(0.25000E+01)*math.sqrt(0.14600E+04)), "uniform"],
		["c_RWAI", unif_margin(4.23245e-7*0.01*math.sqrt(0.15E+01)*math.sqrt(0.14600E+04)), "uniform"],
		["c_SM_act", unif_margin(4.23245e-7*0.01*math.sqrt(2.49)*math.sqrt(0.30000E+6)), "uniform"],
		["c_PM", unif_margin(4.23245e-7*0.01*math.sqrt(0.18860E+02)*math.sqrt(0.77400E+06)), "uniform"],
		["c_PM_act", unif_margin(4.23245e-7*0.01*math.sqrt(0.18860E+02)*math.sqrt(0.30000E+6)), "uniform"],
		["c_petal", unif_margin(4.23245e-7*0.01*math.sqrt(0.18860E+02)*math.sqrt(0.9000E+08)), "uniform"],
		["zeta_sunshield", unif_margin(0.005*5), "uniform"],
		["zeta_isolator", unif_margin(0.005*20), "uniform"],
		["zeta_solarpanel", unif_margin(0.005*20), "uniform"],	
		["Ru", unif_margin(3000), "uniform"],
		["fc", unif_margin(30), "uniform"],
		["Tst", unif_margin(20), "uniform"],
		["Srg", unif_margin(3e-14), "uniform"],
		["Sst", unif_margin(2), "uniform"],
		["Tgs", unif_margin(0.04), "uniform"],
		["lambda_", unif_margin(1e-6), "uniform"],
		["Ro", unif_margin(0.98), "uniform"],
		["QE", unif_margin(0.8), "uniform"],
		["Mgs", unif_margin(15), "uniform"],
		["fca", unif_margin(0.01), "uniform"],
		["Kc", [0,1], "uniform"],
		["Kcf", unif_margin(2000), "uniform"],
	]
	var_names = [pdef[0] for pdef in param_defs]
	var_bounds = [pdef[1] for pdef in param_defs]
	var_dists = [pdef[2] for pdef in param_defs]

	def model(params):
		###First, process theta, sampling using the new precisions:
		theta = []
		x = deepcopy(problem.x_default)
		xrandom = ["Ru","fc","Tst","Srg","Sst","Tgs","lambda_","Ro","QE","Mgs","fca","Kc","Kcf"]
		for i,name in enumerate(var_names):
			if name in xrandom:
				idx = xrandom.index(name)
				x[idx] == params[i]
			else:
				theta.append(params[i])
		
		###Finally, return QoI
		return problem.H(theta, x, verbose=False)
		
	###Generate some new samples and save
	#I want to save samples as often as possible. Therefore, iterate through N two at a time, corresponding to 2(p+2) model evals at a time
	for _ in range(int(N/2)):
		saltelli_eval_sample(filename, 2, var_names, var_dists, var_bounds, model, doPrint=True)

#sensitivity analysis of HLVA
def vv_SA_jitter_evaluate(problem, do_subset):
	var_names = ["Us","Ud","Qc","I_SMhubt","I_SMhuba","K_yPM","I_xRWA","I_yRWA","I_RWt","I_RWa","I_ISOa","I_ISOt","K_yISO","K_xISO","I_bus","I_propt","I_propa","I_i1","I_i2","I_i3","A_sptop","D_sp","t_sp","I_ss","K_rad1","K_rad2","K_rISO","K_act1","K_act2","I_iso","K_zpet","K_pm1","K_pm3","K_pm4","K_pm5","K_pm6","K_act_pm2","K_act_pm3","K_act_pm4","K_act_pm5","K_act_pm6","K_xpet","c_RWA","c_RWAI","c_SM_act","c_PM","c_PM_act","c_petal","zeta_sunshield","zeta_isolator","zeta_solarpanel","Ru","fc","Tst","Srg","Sst","Tgs","lambda_","Ro","QE","Mgs","fca","Kc","Kcf"]

	#S, ST, n_eval = saltelli_indices("SA_jitter", var_names, do_subset=0, doPrint=True)
	total_order_convergence_tests(1200, "SA_jitter", var_names, do_subset=do_subset)

def vv_SA_jitter_convergence(problem):
	var_names = ["Us","Ud","Qc","I_SMhubt","I_SMhuba","K_yPM","I_xRWA","I_yRWA","I_RWt","I_RWa","I_ISOa","I_ISOt","K_yISO","K_xISO","I_bus","I_propt","I_propa","I_i1","I_i2","I_i3","A_sptop","D_sp","t_sp","I_ss","K_rad1","K_rad2","K_rISO","K_act1","K_act2","I_iso","K_zpet","K_pm1","K_pm3","K_pm4","K_pm5","K_pm6","K_act_pm2","K_act_pm3","K_act_pm4","K_act_pm5","K_act_pm6","K_xpet","c_RWA","c_RWAI","c_SM_act","c_PM","c_PM_act","c_petal","zeta_sunshield","zeta_isolator","zeta_solarpanel","Ru","fc","Tst","Srg","Sst","Tgs","lambda_","Ro","QE","Mgs","fca","Kc","Kcf"]
	
	###Perform the analysis at a few evaluation points
	list_S = []
	list_ST = []
	list_n = [10,50,100,500,1000,5000,10000,50000]
	for n in list_n:
		#S, ST, n_eval = saltelli_indices("SA_jitter", var_names, do_subset=0, doPrint=True)
		#list_S.append(S)
		#list_ST.append(ST)
		total_order_convergence_tests(800, "SA_jitter", var_names, do_subset=n**2)

"""
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
"""

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Function to run for this vvopt analysis')
	parser.add_argument('-n', type=int, default=0, help='Number of iterations to give to the function')
	parser.add_argument('--filename', metavar='string', default="SA_jitter", help='Base name to five to SA_jitter_sample')
	args = parser.parse_args()
	
	###Problem Definition
	problem = construct_jwst_jitter_problem()
	theta_nominal = problem.theta_nominal
	
	d0s = [0 for _ in problem.d_names]
	d1s = [1 for _ in problem.d_names]
	d2s = [2 for _ in problem.d_names]
	d3s = [3 for _ in problem.d_names]

	###Uncertainty Quantification
	if args.run == "system_nominal":
		vv_system(problem, theta_nominal)
		
	elif args.run == "experiment_nominal":	
		vv_experiment(problem, theta_nominal, d2s)
	
		"""
	elif args.run == "system_model":
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
		
	elif args.run == "cheap_design":
		print("Given the nominal theta:", theta_nominal)
		print("and the cheapest design:", d_min)
		y_cheap = problem.eta(theta_nominal, d_min, err=True)
		
		print("Cheapest y:", y_cheap)
		print("Cost of design:", problem.G(d_min))
	"""
	
	elif args.run == "UA_theta":
		vv_UA_theta(problem, n=args.n)	

	elif args.run == "UP_QoI":
		vv_UP_QoI(problem, 0, n=args.n)
		
	elif args.run == "UP_QoI_samples":
		vv_UP_QoI_samples(base_name="SA_jitter", doPrint=True, do_subset=args.n)
		
	elif args.run == "UP_exp":
		vv_UP_exp(problem, d_historical, theta_nominal, n=args.n)
	
	elif args.run == "SA_jitter_sample":
		vv_SA_jitter_sample(problem, N=args.n, filename=args.filename)

	elif args.run == "SA_jitter_sample":
		vv_SA_jitter_sample(problem, N=args.n)

	elif args.run == "SA_jitter_evaluate":
		vv_SA_jitter_evaluate(problem, do_subset=args.n)
		
	elif args.run == "SA_jitter_convergence":
		vv_SA_jitter_convergence(problem)
	
		"""
	#Still needs massaging...
	#if args.run == "SA_exp":
	#	vv_SA_exp(problem, d_historical)

	###Optimal Bayesian Experimental Design
	elif args.run == "BN_sample":
		rate = 10 if args.filename=="SA_jitter" else int(args.filename)	
		bn_sampling(problem, savefile="BN_samples", N=args.n, buffer_rate=rate, doPrint=True)
	
	elif args.run == "BN_train":
		#Train the BN off of the saved data
		q, _ = bn_load_samples(problem, savefile="BN_samples", doPrint=True, doDiagnostic=True)
		gmm = bn_train_from_file(problem, savefile="BN_samples", do_subset=args.n, doPrint=True)
		
		#Save the GMM to a file
		#filename = "BN_" + str(len(q)) + '.csv'
		filename = "BN_model.csv"
		bn_save_gmm(gmm, gmm_file=filename)
		
	elif args.run == "BN_evaluate":
		#Run the validation test
		#gmm = bn_load_gmm("BN_model.csv")
		#bn_measure_model_mse(problem, gmm, N=args.n, doPrint=True)
		
		#bn_compare_model_covariance(problem, "BN_samples", "BN_model_100000", doPrint=True)
		
		bn_evaluate_model_likelihood(problem, gmmfile="BN_model_1000", datafile="BN_samples", N_val=0, do_subset=1000, doPrint=True)
		
	elif args.run == "BN_convergence":
		#Run the convergence test
		#bn_measure_stability_convergence(problem, , N_val=args.n, doPrint=True)
		bn_measure_likelihood_convergence(problem, "BN_samples", doPrint=True)
	
	elif args.run == "OBED_test":
		#Load the GMM from file
		gmm = bn_load_gmm("BN_model.csv")
	
		#Calculate U for several different designs
		U_varH_gbi_joint(d_historical, problem, gmm, n_mc=args.n, ncomp=0, doPrint=True)
		U_varH_gbi_joint(d_min, problem, gmm, n_mc=args.n, ncomp=0, doPrint=True)
	"""
	
	else:
		print("I dont recognize the command",args.run)
