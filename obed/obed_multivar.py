#This details an obed script that can be called regardless of the problem

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
import dill

sys.path.append('..')
from obed.mcmc import *
from obed.pdf_estimation import *


"""
This function is part of the OBED problem, solving the utility for a given d
This may be an inefficient brute-force approach
This assumes a utility u(d,y,theta) = 1/Var[H(theta|y,d)], minimizing the variance of the HLVA evaluated over the posterior
A higher function must optimize over the outputs of this model to find d*, the d that maximized U_d
"""
#This is like U_brute_varH, except we're reusing loop 1 samples in loop 3, like Huan & Marzouk
#Here, we're using straightforward Metropolis-Hastings to solve the multivariate MCMC problem
#This calculates the probability of meeting the requirement in distribution - used as constraint for the cost-optimization
def U_probreq(d, problem, mcmc_proposal, prop_width, likelihood_kernel, minreq=None, maxreq=None, n_mc=10**4, n_mcmc=10**3, burnin=0, lag=1, doPrint=False):   
	#requirement: do min, max, or both, but not neither
	if minreq==None and maxreq==None:
		print("U_probreq needs a requirement on the QoI in order to determine probability")
		exit()
	
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	pthetas = problem.prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	
	#Then i want to create an estimated pdf of y as a fn of theta
	#I'll generate thetas uniformly across the domain (theta_domains), and sample ys from that
	#likelihood_kernel
	
	U_list = []
	for i,y in enumerate(Y1_list): #MC loop		
		mcmc_trace,_,_ = mcmc_kernel(y, likelihood_kernel, mcmc_proposal, prop_width, problem.prior_rvs, problem.prior_pdf_unnorm, n_mcmc, burnin, lag, doPlot=False, doPrint=False)
	
		#Now, use my posterior samples to calculate H(theta|y,d) samples
		H_theta_posterior = [problem.H(theta) for theta in mcmc_trace]
		
		#compute a specific probability P(H > req)
		#easily calculated with sample probability, just count em up and find the ratio
		num_true = 0
		for h in H_theta_posterior:
			#check if requirement is satisfied, ignoring None's:
			#if (minreq==None and h <= maxreq) or (h >= minreq and maxreq==None) or (h >= minreq and h <= maxreq):
			if minreq==None:
				num_true += int(h <= maxreq)
			elif maxreq==None:
				num_true += int(h >= minreq)
			else:
				num_true += int(h <= maxreq and h >= minreq)
		prob = num_true / len(H_theta_posterior) #is this the best estimator for a probability?
		u = prob

		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u), flush=True)
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list


def U_probreq_multi(d, problem, mcmc_proposal, prop_width, likelihood_kernel, minreq=None, maxreq=None, n_mc=10**4, n_mcmc=10**3, burnin=0, lag=1, doPrint=False):   
	#requirement: do min, max, or both, but not neither
	if minreq==None and maxreq==None:
		print("U_probreq needs a requirement on the QoI in order to determine probability")
		exit()
	
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	pthetas = problem.prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	
	#Then i want to create an estimated pdf of y as a fn of theta
	#I'll generate thetas uniformly across the domain (theta_domains), and sample ys from that
	#likelihood_kernel
	
	U_list = []
	#Define the mcmc loop as a function of one step
	#see _task_mcmc_y	

	#Run the map that calls all of the n_MC MCMC's
	print("U_probreq_multi running on ",mp.cpu_count(),"processors.")
	with mp.Pool() as pool:
		task = partial(_task_mcmc_y, kernel=likelihood_kernel, proposal=mcmc_proposal, width=prop_width, prob=problem, minreq=minreq, maxreq=maxreq, n=n_mcmc, burnin=burnin, lag=lag)
		U_list = pool.map(task, Y1_list)
	
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list
	
	
def _task_mcmc_y(y, kernel, proposal, width, prob, minreq, maxreq, n, burnin, lag):
	mcmc_trace,_,_ = mcmc_kernel(y, kernel, proposal, width, prob.prior_rvs, prob.prior_pdf_unnorm, n, burnin, lag, doPlot=False, doPrint=False)

	#Now, use my posterior samples to calculate H(theta|y,d) samples
	H_theta_posterior = [prob.H(theta) for theta in mcmc_trace]
	
	#compute a specific probability P(H > req)
	#easily calculated with sample probability, just count em up and find the ratio
	num_true = 0
	for h in H_theta_posterior:
		#check if requirement is satisfied, ignoring None's:
		#if (minreq==None and h <= maxreq) or (h >= minreq and maxreq==None) or (h >= minreq and h <= maxreq):
		if minreq==None:
			num_true += int(h <= maxreq)
		elif maxreq==None:
			num_true += int(h >= minreq)
		else:
			num_true += int(h <= maxreq and h >= minreq)
	prob = num_true / len(H_theta_posterior) #is this the best estimator for a probability?
	u = prob
	print("Ran mcmc on",str(y),"for u =",str(u), flush=True)
	return u
	
def U_probreq_1step(d, problem, mcmc_proposal, prop_width, kernel_pickle, minreq=None, maxreq=None, n_mcmc=10**3, burnin=0, lag=1, doPrint=True):   
	#requirement: do min, max, or both, but not neither
	if minreq==None and maxreq==None:
		print("U_probreq needs a requirement on the QoI in order to determine probability")
		exit()
	
	#Generate one y sampled from likelihood fn, p(y|theta,d)p(theta)
	ptheta = problem.prior_rvs(1)
	Y1 = problem.eta(ptheta, d)
	
	#Then i want to create an estimated pdf of y as a fn of theta
	#I'll generate thetas uniformly across the domain (theta_domains), and sample ys from that
	#unpickle the likelihood kernel:
	with open(kernel_pickle, 'rb') as f:
		likelihood_kernel = dill.load(f)
	
	mcmc_trace,_,_ = mcmc_kernel(Y1, likelihood_kernel, mcmc_proposal, prop_width, problem.prior_rvs, problem.prior_pdf_unnorm, n_mcmc, burnin, lag, doPlot=False, doPrint=False)

	#Now, use my posterior samples to calculate H(theta|y,d) samples
	H_theta_posterior = [problem.H(theta) for theta in mcmc_trace]
	
	#compute a specific probability P(H > req)
	#easily calculated with sample probability, just count em up and find the ratio
	num_true = 0
	for h in H_theta_posterior:
		#check if requirement is satisfied, ignoring None's:
		#if (minreq==None and h <= maxreq) or (h >= minreq and maxreq==None) or (h >= minreq and h <= maxreq):
		if minreq==None:
			num_true += int(h <= maxreq)
		elif maxreq==None:
			num_true += int(h >= minreq)
		else:
			num_true += int(h <= maxreq and h >= minreq)
	prob = num_true / len(H_theta_posterior)
	u = prob

	if doPrint:
		print(str(u)+",",str(Y1),flush=True)
		
	#compute an in-distribution probability
	return u
	
	
def U_probreq_1step_nokernel(d, problem, mcmc_proposal, prop_width, minreq=None, maxreq=None, n_mcmc=10**3, n_pde=1000, burnin=0, lag=1, doPrint=True):   
	#requirement: do min, max, or both, but not neither
	if minreq==None and maxreq==None:
		print("U_probreq needs a requirement on the QoI in order to determine probability")
		exit()
	
	#Generate one y sampled from likelihood fn, p(y|theta,d)p(theta)
	ptheta = problem.prior_rvs(1)
	Y1 = problem.eta(ptheta, d)
	
	mcmc_trace,_,_ = mcmc_multigauss_likelihood(Y1, d, mcmc_proposal, prop_width, problem.eta, problem.prior_rvs, problem.prior_pdf_unnorm, n_mcmc, n_pde, burnin, lag, doPlot=False, doPrint=False)

	#Now, use my posterior samples to calculate H(theta|y,d) samples
	H_theta_posterior = [problem.H(theta) for theta in mcmc_trace]
	
	#compute a specific probability P(H > req)
	#easily calculated with sample probability, just count em up and find the ratio
	num_true = 0
	for h in H_theta_posterior:
		#check if requirement is satisfied, ignoring None's:
		#if (minreq==None and h <= maxreq) or (h >= minreq and maxreq==None) or (h >= minreq and h <= maxreq):
		if minreq==None:
			num_true += int(h <= maxreq)
		elif maxreq==None:
			num_true += int(h >= minreq)
		else:
			num_true += int(h <= maxreq and h >= minreq)
	prob = num_true / len(H_theta_posterior)
	u = prob

	if doPrint:
		print(str(u)+",",str(Y1)+",",str(np.mean(H_theta_posterior))+",",str(np.std(H_theta_posterior, ddof=1)), flush=True)
		
	#compute an in-distribution probability
	return u
