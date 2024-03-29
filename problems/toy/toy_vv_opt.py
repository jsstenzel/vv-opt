import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import csv
import dill

sys.path.append('../..')
#focal plane
from problems.toy.toy_problem import *
#analysis
from obed.obed_multivar import *
from obed.mcmc import *
from obed.pdf_estimation import *
from uq.uncertainty_propagation import *
from uq.sensitivity_analysis import *


################################
#Useful definitions
################################

def proposal_fn_gamma(theta_curr, proposal_width):
	theta_prop = [0] * len(theta_curr)
	for i,_ in enumerate(theta_prop):
		#proposal dists are gammas, to match
		mean = abs(theta_curr[i])
		stddev = proposal_width
		variance = [w**2 for w in stddev]
		alpha = mean**2 / variance
		beta = mean / variance
		theta_prop[i] = scipy.stats.gamma.rvs(size=1, a=alpha, scale=1.0/beta)[0]
	return theta_prop
	
def proposal_fn_norm(theta_curr, proposal_width):
	theta_prop = [0] * len(theta_curr)
	for i,_ in enumerate(theta_prop):
		#proposal dists are gammas, to match
		mean = abs(theta_curr[i])
		stddev = proposal_width[i]
		theta_prop[i] = scipy.stats.norm.rvs(size=1, loc=mean, scale=stddev)[0]
	return theta_prop

toy_minreq = 2 #want H to be at least 2

theta_nominal = [0.5]
QoI_nominal = toy.H(theta_nominal)

d_historical = [.2]
d_worst = [.05]
d_best = [.5]
		
y_nominal = toy_eta(dict(zip(toy.theta_names, theta_nominal)), dict(zip(toy.d_names, d_historical)), dict(zip(toy.x_names, toy.x_default)), err=False)
#print(y_nominal)

################################
#Analysis functions
################################

def toy_vv_nominal():
	print("QoI requirement:", toy_minreq)
	print("Nominal QoI:", QoI_nominal)
	print("")
	print("Historical experimental design d:",d_historical)
	print("Nominal experimental result y:",y_nominal)


###uncertainty analysis
def toy_vv_UP_QoI():
	#uncertainty propagation of HLVA
	uq_thetas = toy.prior_rvs(10000)
	Qs = [toy.H(theta) for theta in uq_thetas]
	uncertainty_prop_plots(uq_thetas, xlabs=toy.theta_names)
	uncertainty_prop_plot(Qs, xlab="QoI")

	#prob of meeting req along priors:
	count_meetreq = 0
	for Q in Qs:
		if Q <= toy_minreq:
			count_meetreq += 1
	prob_meetreq = count_meetreq / len(Qs)
	print("Probability of meeting requirement given priors:", prob_meetreq)

#sensitivity analysis of HLVA
def toy_vv_SA_QoI():
	#it'll be straightforward to see the dependence of QoI on theta
	Si = sobol_saltelli(toy.H, 
						2**5, #SALib wants powers of 2 for convergence
						var_names=toy.theta_names, 
						var_dists=[prior[0] for prior in toy.priors], 
						var_bounds=[prior[1] for prior in toy.priors], 
						conf = 0.95, doSijCalc=False, doPlot=True, doPrint=True)

#Uncertainty analysis of the experiment models
def toy_vv_UP_exp(dd, savefig=False):
	print("Likelihood distribution for nominal historical case:",flush=True)
	tt = toy.prior_rvs(1); print(tt)
	ysample_nominal = [toy.eta(tt, dd) for _ in range(10000)]
	uncertainty_prop_plots(ysample_nominal, xlabs=toy.y_names, saveFig='')
	
	#Also plot the y's generated by the joint distribution p(y|theta,d)p(theta)
	print("Experiments simulated from the joint distribution:",flush=True)
	uq_thetas = toy.prior_rvs(10000)
	uq_ys = [toy.eta(theta, dd) for theta in uq_thetas]
	uncertainty_prop_plots(uq_ys, xlabs=toy.y_names, c='orchid', saveFig='')

"""
#sensitivity analysis of the experiment models
def toy_vv_SA_exp(dd):
		#we want to see sensitivity of each yi to their inputs, theta and x
		#instead of trying to break eta into its constitutent experiments, i think i want to just use all of eta,
		#and doing separate sensitivity analysis on each of the y[i] coming from theta, over entire x and theta span
		#sobol_saltelli expects a function that takes a single list of parameters
		def eta_wrap(param): #gain
			gain = param[0]
			rn = param[1]
			dc = param[2]
			_x = dict(zip(toy.x_names, toy.x_default)) 
			_x["sigma_dc"] = param[3]
			_x["P_signal"] = param[4]
			_x["P_noise"] = param[5]
			_x["T_ccd"] = param[6]
			_x["sigma_E"] = param[7]
			_x["w"] = param[8]
			_x["activity_cd109"] = param[9]
			y = toy_eta(theta, d, x, err=True)
			return y

		expvar_names = toy.theta_names + toy.x_names
		expvar_dists = [prior[0] for prior in toy.priors] + [prior[0] for prior in toy.x_dists]
		expvar_bounds=[prior[1] for prior in toy.priors] + [prior[1] for prior in toy.x_dists]
			
		#so ugly!!!
		Si_1 = sobol_saltelli(eta_wrap, 
							2**8, #SALib wants powers of 2 for convergence
							var_names=expvar_names, 
							var_dists=expvar_dists,
							var_bounds=expvar_bounds,
							conf = 0.99, doSijCalc=False, doPlot=True, doPrint=True)
"""

""" fix for toy
def toy_vv_mcmc_convergence(): #convergence study with mcmc_multivar
	print("Generating kernel",flush=True)
	n_kde_axis = 47 #47^3 equals about 10^5
	kde_gains = np.linspace(0,3,n_kde_axis)
	kde_rn = np.linspace(1,4,n_kde_axis)
	kde_dc = np.linspace(0,.01,n_kde_axis)
	kde_thetas = np.vstack((np.meshgrid(kde_gains, kde_rn, kde_dc))).reshape(3,-1).T #thanks https://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays
	#uncertainty_prop_plots(kde_thetas, xlabs=['gain','rn','dc'])
	kde_ys = [toy.eta(theta, d_historical) for theta in kde_thetas]
	#uncertainty_prop_plots(kde_ys, xlabs=['Y0','Y1','Y2'])
	likelihood_kernel, kde_ythetas = general_likelihood_kernel(kde_thetas, kde_ys)
	kde_plot(likelihood_kernel, kde_ythetas, plotStyle='together', ynames=['gain','rn','dc','y1','y2','y3']) #this guy is 6d? or 3d? dont know how to plot it
	
	print("mcmc convergence study",flush=True)
	#crucially, this is being evaluated at y=y_nominal. Good and representative, but we'll have to check robustness after.
	n_mcmc = 10**5
	prop_fn = proposal_fn_norm
	#prop_width = [.1,.05,.0001] #rough guess from looking at prior + a few short runs - did 10**5 run
	#prop_width = [5.17612361e-05, 2.53953018e-02, 3.24906620e-07] #from covariance of the above study #this led to 89% acceptance rate, too high! #whoops, used variance instead of stddev
	prop_width = [0.0071945282048228735, 0.15935903430820628, 0.0005700058073427601] #stddev from first study
	mcmc_trace, arate, rrate = mcmc_kernel(y_nominal, likelihood_kernel, prop_fn, prop_width, toy.prior_rvs, toy.prior_pdf_unnorm, n_mcmc, burnin=0, lag=1, doPlot=True, legend=toy.theta_names, doPrint=True)
	print(arate, rrate)
	
	#save data, do analysis and plots
	with open('mcmc.csv', 'w', newline='') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter=' ')
		for theta in mcmc_trace:
			csvwriter.writerow(theta)
	print("mean, stddev, covariance of posterior sample:")
	means, stddevs, cov = mcmc_analyze(mcmc_trace,doPlot=True)
	print(means)
	print(stddevs)
	print(cov)
	uncertainty_prop_plot([sample[0] for sample in mcmc_trace], c='limegreen', xlab="Gain [ADU/e-]")
	uncertainty_prop_plot([sample[1] for sample in mcmc_trace], c='limegreen', xlab="Read noise [e-]")
	uncertainty_prop_plot([sample[2] for sample in mcmc_trace], c='limegreen', xlab="Dark current [e-/s]")

#if False: #play with mcmc_nolikelihood
#	#This is really too slow, so I want to stay away from it
#	print("mcmc",flush=True)
#	n_kde = 10**4
#	n_mcmc = 10**2
#	prop_fn = proposal_fn_gamma
#	prop_width = [.3,.3,.002]
#	mcmc_trace, arate, rrate = mcmc_nolikelihood(y_nominal, d_historical, toy.eta, proposal_fn_norm, toy.prior_rvs, toy.prior_pdf_unnorm, n_mcmc, n_kde, burnin=0, lag=1, doPlot=True, legend=toy.theta_names)
#	print("acceptance rate",arate, "randomwalk rate", rrate)
"""
	
""" fix for toy
###mcmc robustness study
def toy_vv_mcmc_robustness():
	#this value should be a good one derived from the above
	prop_width = [0.007167594573520732, 0.17849464019335232, 0.0006344271319903282]

	print("Generating kernel",flush=True)
	n_kde_axis = 47 #47^3 equals about 10^5
	kde_gains = np.linspace(0,3,n_kde_axis)
	kde_rn = np.linspace(1,4,n_kde_axis)
	kde_dc = np.linspace(0,.01,n_kde_axis)
	kde_thetas = np.vstack((np.meshgrid(kde_gains, kde_rn, kde_dc))).reshape(3,-1).T #thanks https://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays
	kde_ys = [toy.eta(theta, d_historical) for theta in kde_thetas]
	likelihood_kernel, kde_ythetas = general_likelihood_kernel(kde_thetas, kde_ys)
	
	print("mcmc convergence robustness study",flush=True)
	for ysample in [toy.eta(tt, d_historical) for tt in toy.prior_rvs(5)]:
		print("ysample mcmc test:",ysample)
		n_mcmc = 6000
		prop_fn = proposal_fn_norm
		mcmc_trace, arate, rrate = mcmc_kernel(ysample, likelihood_kernel, proposal_fn_norm, prop_width, toy.prior_rvs, toy.prior_pdf_unnorm, n_mcmc, burnin=300, lag=1, doPlot=True, legend=toy.theta_names, doPrint=True)
		print(arate, rrate, "			")
		means, stddevs, cov = mcmc_analyze(mcmc_trace,doPlot=True)
		print("mean of posterior sample", means)
		print("stddev of posterior sample", stddevs)
		print("covariance of posterior sample")
		print(cov)
		uncertainty_prop_plot([sample[0] for sample in mcmc_trace], c='limegreen', xlab="Gain [ADU/e-]")
		uncertainty_prop_plot([sample[1] for sample in mcmc_trace], c='limegreen', xlab="Read noise [e-]")
		uncertainty_prop_plot([sample[2] for sample in mcmc_trace], c='limegreen', xlab="Dark current [e-/s]")
"""

""" fix for toy
###pickle the kernel!
def toy_vv_kernel_pickle(dd, name='likelihood_kernel.pkl'):
	print("Generating kernel",flush=True)
	n_kde_axis = 50 #47^3 equals about 10^5
	kde_gains = np.linspace(0,3,n_kde_axis)
	kde_rn = np.linspace(1,4,n_kde_axis)
	kde_dc = np.linspace(0,.01,n_kde_axis)
	kde_thetas = np.vstack((np.meshgrid(kde_gains, kde_rn, kde_dc))).reshape(3,-1).T #thanks https://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays
	kde_ys = [toy.eta(theta, dd) for theta in kde_thetas]
	likelihood_kernel, kde_ythetas = general_likelihood_kernel(kde_thetas, kde_ys)
	
	print("pickle it!")
	with open(name, 'wb') as f:
		dill.dump(likelihood_kernel, f)
"""

""" fix for toy
###obed analysis -- naive multiprocessing approach
def toy_vv_obed_1machine():
	print("Generating kernel",flush=True)
	n_kde_axis = 47 #47^3 equals about 10^5
	kde_gains = np.linspace(0,3,n_kde_axis)
	kde_rn = np.linspace(1,4,n_kde_axis)
	kde_dc = np.linspace(0,.01,n_kde_axis)
	kde_thetas = np.vstack((np.meshgrid(kde_gains, kde_rn, kde_dc))).reshape(3,-1).T #thanks https://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays
	kde_ys = [toy.eta(theta, d_historical) for theta in kde_thetas]
	likelihood_kernel, kde_ythetas = general_likelihood_kernel(kde_thetas, kde_ys)

	print("starting obed",flush=True)
	prop_width = [0.007167594573520732, 0.17849464019335232, 0.0006344271319903282] #stddev from last study

	U, U_list = U_probreq_multi(d_historical, toy, proposal_fn_norm, prop_width, likelihood_kernel, maxreq=3.0, n_mc=1000, n_mcmc=2000, burnin=300, lag=1, doPrint=True)
	print(U)
	print(U_list)
	uncertainty_prop_plot(U_list, c='royalblue', xlab="specific U", saveFig='OBEDresult')
"""	

""" fix for toy
###obed analysis -- robust for sbatch, one step of the monte carlo
def toy_vv_obed_cluster(d, kernel_pkl):
	print("starting obed",flush=True)
	prop_width = [0.007167594573520732, 0.17849464019335232, 0.0006344271319903282] #stddev from last study

	U = U_probreq_1step(d_historical, toy, proposal_fn_norm, prop_width, kernel_pkl, maxreq=3.0, n_mcmc=2000, burnin=300, lag=1, doPrint=True)
"""

def toy_vv_plot_obed_results(file):
	u = []
	H_mean = []
	H_stddev = []
	with open(file) as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			u.append(float(row[0]))
			H_mean.append(float(row[4]))
			H_stddev.append(float(row[5]))
	
	uncertainty_prop_plot(u, c='royalblue', xlab="specific U from MC")
	uncertainty_prop_plot(H_mean, c='royalblue', xlab="H means from MC")
	uncertainty_prop_plot(H_stddev, c='royalblue', xlab="H stddevs from MC")
	print("MC-mean probability of meeting requirement:",np.mean(u))
	print("MC-mean of H_posterior mean:",np.mean(H_mean))
	print("   Probability of MC-mean meeting requirement:",len([1 for hmean in H_mean if hmean>=toy_minreq])/len(H_mean))
	print("MC-mean of H_posterior stddev:",np.mean(H_stddev))
	
	plt.scatter(H_mean, H_stddev)
	plt.show()
	

def toy_vv_test_mcmc_multigauss(yy, dd):
	print("likelihood for ynominal given d_historical, theta_nominal:",eta_multigaussian_logpdf(yy, theta_nominal, d_historical, toy.eta, n_pde=10000))
	print("likelihood for ynominal given d_best, theta_nominal:",eta_multigaussian_logpdf(yy, theta_nominal, d_best, toy.eta, n_pde=10000))
	print("likelihood for ynominal given d_worst, theta_nominal:",eta_multigaussian_logpdf(yy, theta_nominal, d_worst, toy.eta, n_pde=10000))
	
	#prop_width = [.5]
	prop_width = [0.16947063365721085]
	mcmc_trace,acceptance_rate,_ = mcmc_multigauss_likelihood(yy, dd, proposal_fn_norm, prop_width, toy.eta, toy.prior_rvs, toy.prior_pdf_unnorm, n_mcmc=20, n_pde=1000, burnin=0, lag=1, doPlot=True, legend=toy.theta_names, doPrint=True)
	print("Acceptance rate:",acceptance_rate)
	
	print("mean, stddev, covariance of posterior sample:")
	means, stddevs, cov = mcmc_analyze(mcmc_trace,doPlot=True)
	print(means)
	print(stddevs)
	print(cov)
	uncertainty_prop_plots(mcmc_trace, c='limegreen', xlabs=[name+" (posterior)" for name in toy.theta_names])



def toy_vv_obed_nokernel_cluster(dd=d_historical):
	#print("starting obed",flush=True)
	prop_width = [0.17449416794749348] #stddev from multigauss mcmc study

	U = U_probreq_1step_nokernel(d_historical, toy, proposal_fn_norm, prop_width, minreq=toy_minreq, n_mcmc=2000, n_pde=1000, burnin=200, lag=1, doPrint=True)


if __name__ == '__main__':  
	toy_vv_nominal()
	
	toy_vv_UP_QoI()
	
	toy_vv_SA_QoI()
	
	toy_vv_UP_exp(d_historical, False)
	
	toy_vv_SA_exp(d_historical)
	
	"""
	#toy_vv_mcmc_convergence()
	
	#toy_vv_mcmc_robustness()
	
	#toy_vv_kernel_pickle(d_historical, name='likelihood_kernel_dhistorical.pkl')
	#toy_vv_kernel_pickle(d_best, name='likelihood_kernel_dbest.pkl')
	#toy_vv_kernel_pickle(d_worst, name='likelihood_kernel_dworst.pkl')
	
	#toy_vv_obed_1machine()
	
	#toy_vv_obed_cluster(d_historical, "likelihood_kernel.pkl")
	"""
	
	#toy_vv_test_mcmc_multigauss(y_nominal, d_historical)
	
	toy_vv_obed_nokernel_cluster(d_historical)
	
	#print("toy_vv_obed_nokernel_cluster running on ",mp.cpu_count(),"processors.")
	#with mp.Pool() as pool:
	#	do_list = range(10)
	#	OBED_list = pool.map(toy_vv_obed_nokernel_cluster, do_list)

	#toy_vv_plot_obed_results('output_multigauss_fulldata.txt')