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
	
def kde_train(problem, grid_density, d, doDiagnostic=False):
	print("Generating kernel",flush=True)
	kde_gains = scipy.stats.gamma.rvs(size=grid_density, a=1.0999913287843317**2 / 1.5998220769876884e-06**2, scale=1.5998220769876884e-06**2/1.0999913287843317)
	kde_rn = np.linspace(1,4,grid_density)
	kde_dc = np.linspace(0,.01,grid_density)
	kde_thetas = np.vstack((np.meshgrid(kde_rn, kde_dc))).reshape(2,-1).T #thanks https://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays
	if doDiagnostic:
		uncertainty_prop_plots(kde_thetas, xlabs=['rn','dc'])
	kde_ys = [problem.eta(theta, d, x=problem.sample_x(1)) for theta in kde_thetas]
	if doDiagnostic:
		uncertainty_prop_plots(kde_ys, xlabs=['Y1','Y2'])
	likelihood_kernel, kde_ythetas = general_likelihood_kernel(kde_thetas, kde_ys)
	
	print("Pickling kernel",flush=True)
	with open('likelihood_kde.pkl', 'wb') as file:
		pickle.dump(likelihood_kernel, file)
		
	if doDiagnostic:
		kde_plot(likelihood_kernel, kde_ythetas, plotStyle='together', ynames=['rn','dc','y2','y3'])
	
	return likelihood_kernel
	
def fp_prior_update(fp, ydata, d, n_mcmc, loadKDE=False, loadMCMC=False, doDiagnostic=False):
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
		with open('likelihood_kde_seq.pkl', 'rb') as file:
			likelihood_kernel = pickle.load(file)
	else:
		likelihood_kernel = kde_train(fp, grid_density=500, d=d, doDiagnostic=doDiagnostic)
	
	if loadMCMC:
		mcmc_trace = []
		with open('mcmc.csv', 'r', newline='') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=' ')
			for theta in csvreader:
				mcmc_trace.append([float(t) for t in theta])
	else:
		###Then, do MCMC to get posterior samples
		prop_cov = [[1.13951239e-04, 1.67431766e-07],
					[1.67431766e-07, 1.03994093e-07]]
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
		with open('mcmc_seq.csv', 'w', newline='') as csvfile:
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
	
	uncertainty_prop_plot([sample[0] for sample in mcmc_trace], c='limegreen', xlab="Read noise [e-]")
	uncertainty_prop_plot([sample[1] for sample in mcmc_trace], c='limegreen', xlab="Dark current [e-/s]")
	###Lastly, use those samples to plot the posterior predictive distribution
	uncertainty_prop_plot(H_posterior, c='darkgreen', xlab="Posterior QoI: Avg. Noise [e-]", vline=[req])
	#could fit a gaussian to that guy to estimate the probability change!

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
		
	elif args.run == "OBED_test":
		#Load the GMM and presampled y from file
		print("Loading GMM and presamples...",flush=True)
		gmm = bn_load_gmm("BN_sequential_model_4000000_ncomp45.pkl")
		presampled_ylist = bn_load_y(problem, "BN_sequential_samples.csv", doPrint=False, doDiagnostic=False)
		
		#Calculate U for several different designs
		dstar = [14, 21, 4.35234, 2.75273]
		U_dstar, _ = U_varH_gbi_joint_presampled(dstar, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=True)
		print("U_dstar:",U_dstar,flush=True)
		dstarstar = [15,21,14.1764,2.76208]
		U_dstarstar, _ = U_varH_gbi_joint_presampled(dstarstar, problem, gmm, presampled_ylist, n_mc=args.n, doPrint=True)
		print("U_dstarstar:",U_dstarstar,flush=True)
		
	#Find the highest utility design, subject to a cost cap
	elif args.run == "OPT_costcap":
		print(35482.15)
		minimize_with_penalty(
		problem, 
		costcap=problem.G(d_historical), 
		gmm_file="BN_sequential_model_4000000_ncomp45.pkl", 
		ylist_file="BN_sequential_samples.csv",
		n_mc=50000, 
		n_tries=1, 
		x0=[15,22,1.0184,2.76765],
		ftol=1e-8, #0.0006407042183632374,
		penalty=10
		)
	
		"""
   Cost    Utility    n_meas_rn    d_num    d_max    d_pow
-------  ---------  -----------  -------  -------  -------
25830.7    0.06744           14       21  44.7397  2.73757"""
	elif args.run == "prior_update":
		###for the optimal design, d1, assume data that is nominal
		d2 = [14,21,44.7397,2.73757]
		#ydata = [yi*1.1 for yi in y_nominal] #test data, just choosing things on the large end of y
		ydata = problem.eta(theta_nominal, d2, err=False)		

		###see what the benefit would be in the batch timeline, i.e. apply the design all at once
		fp_prior_update(fp2, ydata, d_historical, n_mcmc=args.n, loadKDE=False, doDiagnostic=True)

	else:
		print("I don't recognize the command",args.run)
