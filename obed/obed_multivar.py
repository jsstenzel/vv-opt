#This details an obed script that can be called regardless of the problem

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns

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

