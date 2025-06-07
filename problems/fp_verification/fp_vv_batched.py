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
	
def replot_nsga(design_pts, util_err):	
	pareto = [
		[3021.3,0.0235377,0.100038,1,3,2,1.00001,2.76496],
		[3046.8,0.0229329,0.101129,1,5,2,1.00104,2.88654],
		[3100.57,0.0222279,39.7335,1,3,2,1.00002,0.46484],
		[3286.9,0.0215993,122.701,1,4,2,1.00038,2.47238],
		[3332.4,0.0205716,0.100038,62,2,2,1.00121,0.528267],
		[3342.6,0.0204546,0.100056,64,1,2,1.00121,0.528267],
		[3353.83,0.0196779,0.100056,64,1,4,1.00121,0.528267],
		[3358.93,0.0189023,0.100055,65,2,4,1.00024,0.528267],
		[3399.83,0.018136,0.101077,73,1,4,1.01352,0.551472],
		[3415.08,0.0174439,0.100063,76,1,4,1.01356,0.435721],
		[3453.04,0.0167077,0.108727,80,1,7,1.00086,0.518346],
		[3519.29,0.0158149,0.108727,80,13,7,1.00086,0.578323],
		[3550.71,0.014904,0.109703,91,2,15,1.00001,2.68692],
		[3554.93,0.0142274,0.100112,92,2,15,1.00001,2.68692],
		[3565.62,0.0140648,0.100112,92,2,17,1.00041,2.68692],
		[3724.85,0.0132064,0.109725,84,38,18,1.03787,2.79654],
		[3730.8,0.0126392,0.116272,81,42,18,1.09693,2.77286],
		[3798.5,0.0118342,0.118226,89,43,22,1.09219,2.72553],
		[3879.48,0.0115567,0.475658,96,44,23,1.09403,2.72553],
		[6154.89,0.0114627,34.8071,74,10,21,1.00069,0.528237],
		[6307.84,0.0111663,38.2868,71,14,21,1.08007,0.431858],
		[6396.94,0.0103785,34.8071,80,10,21,1.25751,0.529687],
		[6627.34,0.0095346,38.4666,78,12,23,1.15142,0.431858],
		[6950.76,0.00882877,41.9568,80,5,21,1.22864,0.420593],
		[7050.37,0.00861584,44.1265,78,10,21,1.00167,0.528107],
		[7174.39,0.00790691,41.6345,84,14,21,2.56232,0.423237],
		[7406.25,0.0077057,49.8757,75,23,21,1.00084,0.519905],
		[7570.91,0.0070139,49.8757,78,23,21,1.00274,0.519517],
		[8012.48,0.00678689,54.4345,80,15,21,1.26173,0.497562],
		[8256.26,0.00621329,57.5071,80,14,21,1.26173,0.497562],
		[8938.3,0.00598604,62.3074,84,14,22,2.56232,0.42297],
		[12894.5,0.00597515,101.384,91,1,20,1.00215,2.77225],
		[13806.1,0.00583011,111.293,91,1,20,1.00009,2.7701],
		[14626.5,0.0056716,120.206,91,1,20,1.02236,2.69656],
		[15450.6,0.00542667,127.723,92,1,20,1.00126,2.64808],
		[16263.9,0.00535724,134.962,93,1,20,1.01814,2.64659],
		[16264.4,0.00529055,134.969,93,1,20,1.01814,2.77167],
		[17014.1,0.0051372,144.471,92,1,21,1.09697,2.70443],
		[17738.3,0.00498594,152.263,92,1,21,1.01411,2.69028],
		[18399.9,0.00479681,159.378,92,1,21,1.01366,2.69123],
		[19211.7,0.00471889,168.101,92,1,21,1.09882,2.69047],
		[19976,0.00461837,175.448,92,16,21,1.01021,2.69122],
		[20880.5,0.00459253,188.129,91,1,21,1.00019,2.72462],
		[20962.3,0.00437662,188.129,91,15,22,1.00019,2.76492],
		[22579.2,0.00427538,203.601,92,13,21,1.06234,2.68431],
		[23413.1,0.00416489,212.407,92,15,22,1.00038,2.87492],
		[23886.9,0.00404214,217.498,92,15,22,1.03951,2.79583],
		[24844,0.00402226,230.261,91,15,23,1.03551,2.7738],
		[25433.6,0.00401443,236.73,91,15,22,1.0184,2.76765],
		[25917.1,0.00396416,236.73,93,15,22,1.0184,2.76765]
	]

	plot_nsga2([p[0] for p in pareto], [p[1] for p in pareto], design_pts=design_pts, util_err=util_err, showPlot=True, savePlot=False, logPlotXY=[False,False])

def vv_OPT(problem, gmm_file, ysamples_file, design_pts, tolDelta, util_err, do_hrs, do_min, threads, popSize, nMC, displayFreq, samples):
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
		tolDelta=tolDelta,
		nPeriod=5,
		nMonteCarlo=nMC,
		GMM=gmm,
		Ylist=presampled_ylist,
		displayFreq=displayFreq,
		initial_pop=samples
	)
	plot_nsga2(costs, utilities, design_pts, util_err=util_err, showPlot=True, savePlot=False, logPlotXY=[False,False])

def kde_train(grid_density, d, doDiagnostic=False):
	print("Generating kernel",flush=True)
	kde_gains = np.linspace(0,3,grid_density)
	kde_rn = np.linspace(1,4,grid_density)
	kde_dc = np.linspace(0,.01,grid_density)
	kde_thetas = np.vstack((np.meshgrid(kde_gains, kde_rn, kde_dc))).reshape(3,-1).T #thanks https://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays
	if doDiagnostic:
		uncertainty_prop_plots(kde_thetas, xlabs=['gain','rn','dc'])
	kde_ys = [fp.eta(theta, d_historical) for theta in kde_thetas]
	if doDiagnostic:
		uncertainty_prop_plots(kde_ys, xlabs=['Y0','Y1','Y2'])
	likelihood_kernel, kde_ythetas = general_likelihood_kernel(kde_thetas, kde_ys)
	
	print("Pickling kernel",flush=True)
	with open('likelihood_kde.pkl', 'wb') as file:
		pickle.dump(likelihood_kernel, file)
		
	if doDiagnostic:
		kde_plot(likelihood_kernel, kde_ythetas, plotStyle='together', ynames=['gain','rn','dc','y1','y2','y3'])
	
	return likelihood_kernel

def fp_prior_update(ydata, d, n_mcmc, loadKDE=False, loadMCMC=False, doDiagnostic=False):
	def proposal_fn_norm(theta_curr, prop_cov):
		#theta_prop = [0] * len(theta_curr)
		#for i,_ in enumerate(theta_prop):
		#	#proposal dists are gammas, to match
		#	mean = abs(theta_curr[i])
		#	stddev = proposal_width[i]
		#	theta_prop[i] = scipy.stats.norm.rvs(size=1, loc=mean, scale=stddev)[0]
		
		theta_prop = scipy.stats.multivariate_normal.rvs(mean=theta_curr, cov=prop_cov, size=1)
		return theta_prop
		
	###First, get the KDE
	if loadKDE:
		print("Loading kernel",flush=True)
		with open('likelihood_kde.pkl', 'rb') as file:
			likelihood_kernel = pickle.load(file)
	else:
		likelihood_kernel = kde_train(grid_density=100, d=d, doDiagnostic=doDiagnostic)
	
	if loadMCMC:
		mcmc_trace = []
		with open('mcmc.csv', 'r', newline='') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=' ')
			for theta in csvreader:
				mcmc_trace.append([float(t) for t in theta])
	else:
		###Then, do MCMC to get posterior samples
		prop_cov = [[ 0.2**2, 0, 0],
					[ 0, 0.25**2, 0],
					[ 0, 0, 0.001**2]]
		prop_cov = [[ 9.56717238e-05, -1.74880135e-04,  1.04860536e-07],
					[-1.74880135e-04,  3.81200241e-02, -5.36382479e-06],
					[ 1.04860536e-07, -5.36382479e-06,  5.69947631e-07]]
		mcmc_trace, arate, rrate, last_prop_width = mcmc_kernel(	
			ydata, 
			likelihood_kernel, 
			proposal_fn_norm, 
			prop_cov, 
			fp.prior_rvs, 
			fp.prior_pdf_unnorm, 
			n_mcmc, 
			burnin=300, 
			lag=1, 
			doAdaptive=100,
			doPlot=True, 
			legend=fp.theta_names, 
			doPrint=True)
		print(arate, rrate)
		print(last_prop_width)
		
		###Analyze the mcmc posterior samples and plot the posterior distributions
		#save data, do analysis and plots
		with open('mcmc.csv', 'w', newline='') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter=' ')
			for theta in mcmc_trace:
				csvwriter.writerow(theta)
		
	###Analyze the mcmc posterior samples and plot the posterior distributions
	#save data, do analysis and plots
	print("mean, stddev, covariance of posterior sample:")
	means, stddevs, cov = mcmc_analyze(mcmc_trace,doPlot=True)
	print(means)
	print(stddevs)
	print(cov)
	H_posterior = [fp.H(tt) for tt in mcmc_trace]
	print("QoI posterior mean:", np.mean(H_posterior))
	print("QoI posterior stddev:", np.std(H_posterior))
	print("Posterior probability of meeting the requirement: ", np.sum([int(h <= req) for h in H_posterior])/len(H_posterior))
	
	uncertainty_prop_plot([sample[0] for sample in mcmc_trace], c='limegreen', xlab="Gain [ADU/e-]")
	uncertainty_prop_plot([sample[1] for sample in mcmc_trace], c='limegreen', xlab="Read noise [e-]")
	uncertainty_prop_plot([sample[2] for sample in mcmc_trace], c='limegreen', xlab="Dark current [e-/s]")
	###Lastly, use those samples to plot the posterior predictive distribution
	uncertainty_prop_plot(H_posterior, c='darkgreen', xlab="Posterior QoI: Avg. Noise [e-]", vline=[req])
	#could fit a gaussian to that guy to estimate the probability change!

def gain_prior_update(ygain, d, n_mcmc, n_kde, loadKDE=False, loadMCMC=False, doDiagnostic=False):
	def proposal_fn_norm(theta_curr, prop_width):		
		theta_prop = scipy.stats.norm.rvs(loc=theta_curr, scale=prop_width, size=1)[0]
		return theta_prop
		
	###First, get the KDE
	if loadKDE:
		print("Loading kernel",flush=True)
		with open('gain_likelihood_kde.pkl', 'rb') as file:
			likelihood_kernel = pickle.load(file)
	else:
		print("Generating kernel",flush=True)
		kde_gains = scipy.stats.uniform.rvs(size=n_kde, loc=0, scale=3-0)
		kde_rn = scipy.stats.gamma.rvs(size=n_kde, a=2.5**2 / 0.25**2, scale=0.25**2/2.5)
		kde_dc = scipy.stats.gamma.rvs(size=n_kde, a=0.001**2 / .001**2, scale=.001**2/0.001)
		t_d = d[0]
		I_d = d[1]
		x = problem.x_dict
		kde_thetas = [[gain,rn,dc] for gain,rn,dc in zip(kde_gains,kde_rn,kde_dc)]
		if doDiagnostic:
			uncertainty_prop_plots(kde_thetas, xlabs=['gain','rn','dc'])
		kde_ys = [gain_exp(theta[0], theta[1], theta[2], t_d, I_d, x, err=True) for theta in kde_thetas]
		if doDiagnostic:
			uncertainty_prop_plot(kde_ys, xlab=['Y gain'])
			print(kde_gains)
			print(kde_ys)
		y_theta_values = np.array([[gain,y] for gain,y in zip(kde_gains,kde_ys)])
		likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)
		
		print("Pickling kernel",flush=True)
		with open('gain_likelihood_kde.pkl', 'wb') as file:
			pickle.dump(likelihood_kernel, file)
			
		if doDiagnostic:
			kde_plot(likelihood_kernel, y_theta_values, plotStyle='together', ynames=['gain','y1'])
	
	if loadMCMC:
		mcmc_trace = []
		with open('gain_mcmc.csv', 'r', newline='') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=' ')
			for theta in csvreader:
				mcmc_trace.append([float(t) for t in theta])
	else:
		###Then, do MCMC to get posterior samples
		prop_width = 0.2**2
		def gain_rvs(n):
			return scipy.stats.gamma.rvs(size=1, a=1.1**2 / 0.2**2, scale=0.2**2/1.1)[0]
			
		def gain_pdf(g):
			return scipy.stats.gamma.pdf(g, a=1.1**2 / 0.2**2, scale=0.2**2/1.1)
			
		mcmc_trace, arate, rrate, last_prop_width = mcmc_kernel_1D(	
			ygain, 
			likelihood_kernel, 
			proposal_fn_norm, 
			prop_width, 
			gain_rvs, 
			gain_pdf, 
			n_mcmc, 
			burnin=300, 
			lag=1, 
			doAdaptive=100,
			doPlot=True, 
			legend=["gain"], 
			doPrint=True)
		print(arate, rrate)
		print(last_prop_width)
		
		###Analyze the mcmc posterior samples and plot the posterior distributions
		#save data, do analysis and plots
		with open('gain_mcmc.csv', 'w', newline='') as csvfile:
			csvwriter = csv.writer(csvfile, delimiter=' ')
			for theta in mcmc_trace:
				csvwriter.writerow([theta])
		
	###Analyze the mcmc posterior samples and plot the posterior distributions
	#save data, do analysis and plots
	print("mean, variance of posterior sample:")
	print(np.mean(mcmc_trace))
	print(np.var(mcmc_trace, ddof=1),flush=True)

	uncertainty_prop_plot(mcmc_trace, c='limegreen', xlab="Gain [ADU/e-]")
	###Lastly, use those samples to plot the posterior predictive distribution

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
		1,	  #n_meas_rn
		8,	  #d_num
		9600, #d_max
		2		#d_pow   #approx
	]

	d_max = [
		600,   #t_gain
		100,   #I_gain
		50,	 #n_meas_rn
		100,	#d_num
		12000, #d_max
		3	 #d_pow   #approx
	]	

	d_min = [
		1,   #t_gain
		1,   #I_gain
		1,	  #n_meas_rn
		2,	  #d_num
		1, #d_max
		0.1	  #d_pow   #approx
	]

	theta_nominal = [1.1, 2.5, .001]
	QoI_nominal = fp.H(theta_nominal)
	y_nominal = problem.eta(theta_nominal, d_historical, err=False)
	
	"""
	histcost = fp.G(d_historical)
	bestcost = fp.G([14.34365321,  3, 19, 24, 1.00236047,  2.24834641])
	print("historical design:",histcost, "seconds or",histcost/(3600),"hours")
	print("best design:",bestcost, "seconds or",bestcost/(3600),"hours")
	"""
	design_pts = [
		[problem.G(d_historical),  0.03470218984967855, "d_hist", 0.0003478985402516399],
		[problem.G(d_min), 0.028370595422594045, "d_min", 0.0003478985402516399],
		[problem.G(d_max), 0.09897076229757248, "d_max", 0.0003478985402516399],
	]

	###Uncertainty Quantification
	if args.run == "nominal":
		vv_nominal(problem, req, theta_nominal, y_nominal)
		
	if args.run == "UP_jitter_from_BN":
		Qs, _ = bn_load_samples(problem, savefile="BN_40k_samples", doPrint=True, doDiagnostic=True)
		uncertainty_prop_plot(Qs, xlab="QoI: Avg. Noise [e-]", vline=[req])

		count_meetreq = 0
		for Q in Qs:
			if Q <= req:
				count_meetreq += 1
		prob_meetreq = count_meetreq / len(Qs)
		print("Estimated probability of meeting requirement given priors:", prob_meetreq)
	
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
		U_dmax, _ = U_varH_gbi_joint_presampled(d_max, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=True)
		print("U_dmax:",U_dmax,flush=True)
		U_dmin, _ = U_varH_gbi_joint_presampled(d_min, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=True)
		print("U_dmin:",U_dmin,flush=True)		
	
		
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

	elif args.run == "OPT_test":
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
		second_run = [
		[0, 0, 1, 1, 1, 2, 1, 0.1], #d_min
		[3021.3,0.0235377,0.100038,1,3,2,1.00001,2.76496],
		[3046.8,0.0229329,0.101129,1,5,2,1.00104,2.88654],
		[3100.57,0.0222279,39.7335,1,3,2,1.00002,0.46484],
		[3286.9,0.0215993,122.701,1,4,2,1.00038,2.47238],
		[3332.4,0.0205716,0.100038,62,2,2,1.00121,0.528267],
		[3342.6,0.0204546,0.100056,64,1,2,1.00121,0.528267],
		[3353.83,0.0196779,0.100056,64,1,4,1.00121,0.528267],
		[3358.93,0.0189023,0.100055,65,2,4,1.00024,0.528267],
		[3399.83,0.018136,0.101077,73,1,4,1.01352,0.551472],
		[3415.08,0.0174439,0.100063,76,1,4,1.01356,0.435721],
		[3453.04,0.0167077,0.108727,80,1,7,1.00086,0.518346],
		[3519.29,0.0158149,0.108727,80,13,7,1.00086,0.578323],
		[3550.71,0.014904,0.109703,91,2,15,1.00001,2.68692],
		[3554.93,0.0142274,0.100112,92,2,15,1.00001,2.68692],
		[3565.62,0.0140648,0.100112,92,2,17,1.00041,2.68692],
		[3724.85,0.0132064,0.109725,84,38,18,1.03787,2.79654],
		[3730.8,0.0126392,0.116272,81,42,18,1.09693,2.77286],
		[3798.5,0.0118342,0.118226,89,43,22,1.09219,2.72553],
		[3879.48,0.0115567,0.475658,96,44,23,1.09403,2.72553],
		[6154.89,0.0114627,34.8071,74,10,21,1.00069,0.528237],
		[6307.84,0.0111663,38.2868,71,14,21,1.08007,0.431858],
		[6396.94,0.0103785,34.8071,80,10,21,1.25751,0.529687],
		[6627.34,0.0095346,38.4666,78,12,23,1.15142,0.431858],
		[6950.76,0.00882877,41.9568,80,5,21,1.22864,0.420593],
		[7050.37,0.00861584,44.1265,78,10,21,1.00167,0.528107],
		[7174.39,0.00790691,41.6345,84,14,21,2.56232,0.423237],
		[7406.25,0.0077057,49.8757,75,23,21,1.00084,0.519905],
		[7570.91,0.0070139,49.8757,78,23,21,1.00274,0.519517],
		[8012.48,0.00678689,54.4345,80,15,21,1.26173,0.497562],
		[8256.26,0.00621329,57.5071,80,14,21,1.26173,0.497562],
		[8938.3,0.00598604,62.3074,84,14,22,2.56232,0.42297],
		[12894.5,0.00597515,101.384,91,1,20,1.00215,2.77225],
		[13806.1,0.00583011,111.293,91,1,20,1.00009,2.7701],
		[14626.5,0.0056716,120.206,91,1,20,1.02236,2.69656],
		[15450.6,0.00542667,127.723,92,1,20,1.00126,2.64808],
		[16263.9,0.00535724,134.962,93,1,20,1.01814,2.64659],
		#[16264.4,0.00529055,134.969,93,1,20,1.01814,2.77167],
		[17014.1,0.0051372,144.471,92,1,21,1.09697,2.70443],
		[17738.3,0.00498594,152.263,92,1,21,1.01411,2.69028],
		[18399.9,0.00479681,159.378,92,1,21,1.01366,2.69123],
		[19211.7,0.00471889,168.101,92,1,21,1.09882,2.69047],
		[19976,0.00461837,175.448,92,16,21,1.01021,2.69122],
		[20880.5,0.00459253,188.129,91,1,21,1.00019,2.72462],
		[20962.3,0.00437662,188.129,91,15,22,1.00019,2.76492],
		[22579.2,0.00427538,203.601,92,13,21,1.06234,2.68431],
		[23413.1,0.00416489,212.407,92,15,22,1.00038,2.87492],
		[23886.9,0.00404214,217.498,92,15,22,1.03951,2.79583],
		[24844,0.00402226,230.261,91,15,23,1.03551,2.7738],
		#[25433.6,0.00401443,236.73,91,15,22,1.0184,2.76765],
		[25917.1,0.00396416,236.73,93,15,22,1.0184,2.76765],
		[0, 0, 600, 100, 50, 100, 12000, 3] #d_best
		]

		conf95 = 0.0003478985402516399 #this is 0.5% of the dist between umax and umin
		std_frac = conf95 / (1.96*(0.09829097382501674 - 0.028210891335003707))
		vv_OPT(
			problem,
			gmm_file="BN_batch_model_4000000_ncomp45.pkl",
			ysamples_file="BN_40k_samples.csv",
			design_pts=design_pts,
			tolDelta=std_frac,
			util_err=conf95,
			do_hrs = 0,
			do_min = 0,
			threads = 25,#1 if args.n == 0 else args.n,
			popSize=50,#30 if args.n==0 else args.n,
			nMC=220000,
			displayFreq=5,
			samples=[pop[2:] for pop in second_run],
		)
		
	elif args.run == "OPT_replot":
		replot_nsga(design_pts, 0.0003478985402516399)
		
	#Find the highest utility design, subject to a cost cap
	elif args.run == "OPT_costcap":
		print(problem.G(d_historical))
		minimize_with_penalty(
		problem, 
		costcap=problem.G(d_historical), 
		gmm_file="BN_batch_model_4000000_ncomp45.pkl", 
		ylist_file="BN_40k_samples.csv",
		n_mc=50000, 
		n_tries=1, 
		x0=[236.73,93,15,22,1.0184,2.76765],
		ftol=1e-8, #0.0006407042183632374,
		penalty=10
		)

	elif args.run == "OPT_dumbsearch":
		#Load the GMM and presampled y from file
		print("Loading GMM and presamples...",flush=True)
		gmm = bn_load_gmm("BN_batch_model_4000000_ncomp45.pkl")
		presampled_ylist = bn_load_y(problem, "BN_40k_samples.csv", doPrint=False, doDiagnostic=False)

		while(True):
			di = problem.sample_d(1)
			costi = problem.G(di)
			if 25583.3 < costi and costi <= 35482.15:
				print("Trying a d...",flush=True)
				Ui, _ = U_varH_gbi_joint_presampled(di, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=False)
				if Ui < 0.03470218984967855:
					print(costi, Ui, di, flush=True)
				else:
					print("No good",flush=True)

	elif args.run == "OPT_check":
		dstar_a = [242.052,92,15,22,1.01794,2.73681]
		dstar_b = [235.632,92,14,21,4.35234,2.75273]

		print("Loading GMM and presamples...",flush=True)
		gmm = bn_load_gmm("BN_batch_model_4000000_ncomp45.pkl")
		presampled_ylist = bn_load_y(problem, "BN_40k_samples.csv", doPrint=False, doDiagnostic=False)

		Ua, _ = U_varH_gbi_joint_presampled(dstar_a, problem, gmm, presampled_ylist, n_mc=400000, doPrint=True)
		Ub, _ = U_varH_gbi_joint_presampled(dstar_b, problem, gmm, presampled_ylist, n_mc=400000, doPrint=True)
		print(problem.G(dstar_a), Ua, dstar_a)
		print(problem.G(dstar_b), Ub, dstar_b)

		
	elif args.run == "prior_update":
		###for the optimal design, d1, assume data that is nominal
		d1 = [235.632, 92, 14, 21, 4.35234, 2.75273]
		#ydata = [yi*1.1 for yi in y_nominal] #test data, just choosing things on the large end of y
		ydata = problem.eta(theta_nominal, d1, err=False)		

		###see what the benefit would be in the batch timeline, i.e. apply the design all at once
		fp_prior_update(ydata, d_historical, n_mcmc=args.n, loadKDE=True, loadMCMC=True, doDiagnostic=True)
	
	elif args.run == "sequential_setup":
		d1 = [235.632, 92, 14, 21, 4.35234, 2.75273]
		###conduct gain experiment using d1, assume nominal
		ydata = problem.eta(theta_nominal, d1, err=False)
		gaindata = ydata[0]

		#update only the prior for the gain!
		gain_prior_update(gaindata, d1, n_mcmc=args.n, n_kde=1000000, loadKDE=False, loadMCMC=False, doDiagnostic=True)

	else:
		print("I don't recognize the command",args.run)
