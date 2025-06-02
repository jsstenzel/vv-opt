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
from problems.fp_verification.fp_problem import *
#analysis
from obed.obed_multivar import *
from obed.mcmc import *
from obed.obed_gbi import *
from obed.pdf_estimation import *
from uq.uncertainty_propagation import *
#from uq.sensitivity_analysis import *
from inference.bn_modeling import *
from inference.bn_evaluation import *
from opt.nsga import *



################################
#Analysis functions
################################

def vv_nominal(problem, req, theta_nominal, y_nominal):
	print("QoI requirement:", req)
	QoI_nominal = problem.H(theta_nominal, verbose=True)
	print("Given the nominal theta:", theta_nominal)
	#print("Nominal y:", y_nominal)
	print("Nominal QoI:", QoI_nominal)


def vv_obed_gbi(problem, d):
	U, U_list = U_varH_gbi(d, problem, n_mc=10**5, n_gmm=10**5, doPrint=True)
	print(U)
	#print(U_list)
	uncertainty_prop_plot(U_list, c='royalblue', xlab="specific U")#, saveFig='OBEDresult')
	return U
	
def optimality_study():
	#U_hist = fp_vv_obed_gbi(d_historical)
	U_hist = 0.4981609760099987
	C_hist = np.log(fp.G(d_historical))
	print(C_hist, flush=True)
	
	#U_best = fp_vv_obed_gbi(d_best)
	U_best = 0.48310825880372715
	C_best = np.log(fp.G(d_best))
	
	#U_worst = fp_vv_obed_gbi(d_worst)
	U_worst = 1.2584888138976271
	C_worst = np.log(fp.G(d_worst))
		
	pts = [[C_hist, U_hist, "d_hist"], [C_best, U_best, "d_max"], [C_worst, U_worst, "d_min"]]
	plot_nsga2([C_hist], [U_hist], design_pts=pts, showPlot=True, savePlot=False, logPlotXY=[False,False])
	
	#costs, utilities, designs = ngsa2_problem(fp, nGenerations=200, popSize=30, nMonteCarlo=5*10**3, nGMM=10**3)
	#costs, utilities, designs = ngsa2_problem(fp, nGenerations=10, popSize=2, nMonteCarlo=5*10**3, nGMM=10**3)
	#costs, utilities, designs = ngsa2_problem_parallel(8, fp, hours=0, minutes=0, popSize=12, nMonteCarlo=5*10**3, nGMM=10**3)
	costs, utilities, designs = ngsa2_problem_parallel(8, fp, hours=0, minutes=0, popSize=12, nMonteCarlo=5*10**3, nGMM=5*10**3)

	pts = [[C_hist, U_hist, "d_hist"], [C_best, U_best, "d_max"], [C_worst, U_worst, "d_min"]]
	plot_nsga2(costs, utilities, design_pts=pts, showPlot=True, savePlot=False, logPlotXY=[False,False])
	
def replot_optimal_points():
	#U_hist = fp_vv_obed_gbi(d_historical)
	U_hist = 0.4950101158
	C_hist = np.log(fp.G(d_historical))
	
	#U_best = fp_vv_obed_gbi(d_best)
	U_best = 0.48206003
	C_best = np.log(fp.G(d_best))
	
	#U_worst = fp_vv_obed_gbi(d_worst)
	U_worst = 1.6773693263
	C_worst = np.log(fp.G(d_worst))
	
	#t_gain
	#I_gain
	#n_meas_rn
	#d_num
	#d_max
	#d_pow
	d1 = [34.24242879,  1, 6,  24, 1.00236049,  2.17136203]
	d2 = [14.34365321,  3, 19, 24, 1.00236047,  2.24834641]
	d3 = [9.03012734,   1, 32, 2,  10.6248672,  0.251100910]

	U1 = 0.48891192#fp_vv_obed_gbi(d1)
	C1 = np.log(fp.G(d1))
	U2 = 0.49705924#fp_vv_obed_gbi(d2)
	C2 = np.log(fp.G(d2))
	U3 = 0.53882345#fp_vv_obed_gbi(d3)
	C3 = np.log(fp.G(d3))
	
	pts = [[C_hist, U_hist, "d_hist"], 
			[C_best, U_best, "d_max"], 
			[C_worst, U_worst, "d_min"],
			[C1, U1, "d_1"],
			[C2, U2, "d_2"],
			[C3, U3, "d_3"]
			]
	
	plot_nsga2([C_hist], [U_hist], design_pts=pts, showPlot=True, savePlot=False, logPlotXY=[False,False])

def vv_OPT(problem, gmm_file, ysamples_file, design_pts, epsilon, util_err, do_hrs, do_min, threads, popSize, nMC, displayFreq=10):
	#Load the GMM and presampled y from file
	print("Loading GMM and presamples...",flush=True)
	gmm = bn_load_gmm(gmm_file)
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
		GMM=gmm,
		Ylist=presampled_ylist,
		displayFreq=displayFreq
		#initial_pop=nmcp4_final
	)
	plot_nsga2(costs, utilities, design_pts, util_err=util_err, showPlot=True, savePlot=False, logPlotXY=[False,False])

	

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Functions to run for this vvopt analysis')
	parser.add_argument('-n', type=int, default=0, help='Number of iterations to give to the function')
	parser.add_argument('-ncomp', type=int, default=0, help='Number of components for GMM')	
	args = parser.parse_args()
	
	###Problem Definition
	problem = fp
	req = 4.38 #max noise

	d_historical = [
		20,   #t_gain
		30,   #I_gain
		1,      #n_meas_rn
		8,      #d_num
		9600, #d_max
		2        #d_pow   #approx
	]

	d_best = [
		600,   #t_gain
		100,   #I_gain
		50,     #n_meas_rn
		100,    #d_num
		12000, #d_max
		3     #d_pow   #approx
	]	

	d_worst = [
		1,   #t_gain
		1,   #I_gain
		1,      #n_meas_rn
		2,      #d_num
		1, #d_max
		0.1      #d_pow   #approx
	]

	theta_nominal = [1.1, 2.5, .001]
	QoI_nominal = fp.H(theta_nominal)
	y_nominal = problem.eta(theta_nominal, d_historical, err=False)

	###Uncertainty Quantification
	if args.run == "nominal":
		vv_nominal(problem, req, theta_nominal, y_nominal)
	
	###Train Bayesian network model
	elif args.run == "BN_sample":
		rate = 100
		bn_sampling(problem, savefile="BN_batch_samples", N=args.n, buffer_rate=rate, doPrint=True)
	
	elif args.run == "BN_train":
		#Train the BN off of the saved data
		ncomp = args.ncomp
		q, _ = bn_load_samples(problem, savefile="BN_40k_samples", doPrint=True, doDiagnostic=True)
		gmm = bn_train_from_file(problem, savefile="BN_40k_samples", do_subset=args.n, ncomp=ncomp, doPrint=True)
		
		#Save the GMM to a file
		filename = "BN_batch_model_" + str(len(q)) + "_ncomp" + str(ncomp) + '.pkl'
		#filename = "BN_model.csv"
		bn_save_gmm(gmm, gmm_file=filename)
		
	elif args.run == "BN_examine":
		print("U_dmin:",U_dmin,flush=True)
		bn_compare_model_covariance(problem, "BN_batch_samples", "BN_batch_model_1639027_ncomp200.csv", doPrint=True)
		#bn_plot_data_density(problem, "BN_samples_1639027", "BN_model_1639027_ncomp200.csv", do_subset=100, doGMMPlot=False, doPrint=True)

	elif args.run == "BN_find_ntrain":
		ncomp=20
		N_list = [10**3, 4*10**3, 10**4, 4*10**4, 10**5, 4*10**5]#, 10**6, 4*10**6, 10**7, 4*10**7]
		bn_measure_validation_convergence(problem, "BN_batch_samples", ncomp=ncomp, N_list=N_list, N_val=1000, doPrint=True, doPlot=True)
		#TODO I think i want that fn to produce a plot with many converging lines, but right now it only makes 1 line i think 
	
	elif args.run == "BN_find_ncomp":
		N_list = []
		bn_train_evaluate_ncomp(problem, trainfile="BN_batch_samples", do_subset=4*10**6, doPlot=False, doPrint=True)
		#bn_train_evaluate_ncomp_plot([],[])
		#bn_train_evaluate_ncomp_sanitycheck(problem, trainfile="BN_samples", valfile="BN_validation", doPlot=True, doPrint=True)

	###Optimal Bayesian Experimental Design
	elif args.run == "OBED_test":
		#Load the GMM and presampled y from file
		print("Loading GMM and presamples...",flush=True)
		gmm = bn_load_gmm("BN_batch_model_4000000_ncomp45.pkl")
		presampled_ylist = bn_load_y(problem, "BN_batch_samples.csv", doPrint=False, doDiagnostic=False)
		
		#Calculate U for several different designs
		U_dhist, _ = U_varH_gbi_joint_presampled(d_historical, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=True)
		print("U_dhist:",U_dhist,flush=True)	
		U_dbest, _ = U_varH_gbi_joint_presampled(d_best, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=True)
		print("U_dbest:",U_dbest,flush=True)
		U_dworst, _ = U_varH_gbi_joint_presampled(d_worst, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=True)
		print("U_dworst:",U_dworst,flush=True)		
	
		
	elif args.run == "OBED_convergence":
		#Load the GMM and presampled y from file
		print("Loading GMM and presamples...",flush=True)
		gmm = bn_load_gmm("BN_batch_model_4000000_ncomp45.pkl")
		presampled_ylist = bn_load_y(problem, "BN_batch_samples.csv", doPrint=False, doDiagnostic=False)
		
		#Calculate U_hist for large n_mc, and save the individual MC results
		U_hist,u_1m_list = U_varH_gbi_joint_presampled(d_historical, problem, gmm, presampled_ylist, n_mc=1000000, doPrint=True)
		import pickle
		with open("u_1m_list.pkl", 'wb') as file:
			pickle.dump(u_1m_list, file)
		
	elif args.run == "OBED_convergence_eval":
		import pickle
		with open("u_1m_list.pkl", 'rb') as file:
			u_1m_list = pickle.load(file)
			
		#Take slices of that data for increasing n
		mc_plot_trace_bootstrap(u_1m_list, 60, doLog=False, savePlot=True, doEvery=10000)
	
	"""
	histcost = fp.G(d_historical)
	bestcost = fp.G([14.34365321,  3, 19, 24, 1.00236047,  2.24834641])
	print("historical design:",histcost, "seconds or",histcost/(3600),"hours")
	print("best design:",bestcost, "seconds or",bestcost/(3600),"hours")
	"""
	design_pts = [
		[problem.G(d_historical), 0.034652027998787686, "d_hist", 0.00018466215028162876],
		[problem.G(d_best), 0.09829097382501674, "d_hist", 0.00018466215028162876],
		[problem.G(d_worst), 0.028210891335003707, "d_hist", 0.00018466215028162876],
	]

	
	if args.run == "OBED_plot":
		#U_hist = fp_vv_obed_gbi(d_historical)
		U_hist = 0.4981609760099987
		C_hist = np.log(fp.G(d_historical))
		
		#U_best = fp_vv_obed_gbi(d_best)
		U_best = 0.48310825880372715
		C_best = np.log(fp.G(d_best))
		
		#U_worst = fp_vv_obed_gbi(d_worst)
		U_worst = 1.2584888138976271
		C_worst = np.log(fp.G(d_worst))
		
		pts = [[C_hist, U_hist, "d_hist"], [C_best, U_best, "d_max"], [C_worst, U_worst, "d_min"]]
		plot_ngsa2([C_hist], [U_hist], design_pts=pts, showPlot=True, savePlot=False, logPlotXY=[False,False])

	if args.run == "OPT_test":
		vv_OPT(
			problem, 
			gmm_file="ncomp_testing/BN_model_1639027_ncomp20.pkl", 
			ysamples_file="BN_batch_samples.csv", 
			design_pts=design_pts,
			epsilon=0.01,
			util_err=0.001,
			do_hrs = 0,
			do_min = 0,
			threads = 1 if args.n==0 else args.n,
			popSize=40,# if args.n==0 else args.n,
			nMC=10,
			displayFreq=10
		)

	elif args.run == "OPT":
		conf95 = 0.0003478985402516399
		std_frac = conf95 / (1.96*(0.0047856764435419575 - 0.0022957744137691916))
		vv_OPT(
			problem,
			gmm_file="BN_batch_model_4000000_ncomp45.pkl",
			ysamples_file="BN_batch_samples.csv",
			design_pts=design_pts,
			epsilon=0.001,
			util_err=conf95,
			do_hrs = 0,
			do_min = 0,
			threads = 25,#1 if args.n == 0 else args.n,
			popSize=50,#30 if args.n==0 else args.n,
			nMC=10**4,
			displayFreq=5
		)
