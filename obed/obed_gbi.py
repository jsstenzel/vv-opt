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
from inference.goal_based_inference import *

"""
This function solves a utility U(d,y,theta) = 1/Var[H(theta|y,d)] using goal-based inference
For each d, we train a GMM on the joint distribution of data y and QoI H
For each data point y in the Monte Carlo, we condition the GMM on the data to estimate the posterior H
And we find the variance of that posterior H, wich is our u that we sum up
"""
#This is like U_brute_varH, except we're reusing loop 1 samples in loop 3, like Huan & Marzouk
#Here, we're using straightforward Metropolis-Hastings to solve the multivariate MCMC problem
#This calculates the probability of meeting the requirement in distribution - used as constraint for the cost-optimization
def U_varH_gbi(d, problem, n_mc=10**5, n_gmm=10**4, ncomp=0, doPrint=False):   
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	pthetas = problem.prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	
	#Then we train the GMM offline
	if doPrint:
			print("Training the GMM...", flush=True)
	theta_train = problem.prior_rvs(n_gmm)
	qoi_train = [problem.H(theta) for theta in theta_train]
	y_train = [problem.eta(theta, d) for theta in theta_train]
	gmm = gbi_train_model(theta_train, qoi_train, y_train, ncomp=ncomp)
	
	U_list = []
	for i,y in enumerate(Y1_list): #MC loop		
		betas, mus, Sigs = gbi_condition_model(gmm, y)
	
		#Now, use my posterior predictive to calculate the utility
		H_var = gbi_gmm_variance(betas, mus, Sigs)
		u = H_var
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list
	
def U_varH_gbi_joint(d, problem, gmm_ydq, n_mc=10**5, doPrint=False):   
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	pthetas = problem.prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	print(Y1_list, flush=True)
	
	#We expact we have already trained a joint gmm model p(y,d,q) offline,
	#and will condition it for p(q|y,d)
	
	U_list = []
	for i,y in enumerate(Y1_list): #MC loop		
		vi = np.array(y + d)
	
		#Now, use my posterior predictive to calculate the utility
		H_var = gbi_var_of_conditional_pp(gmm_ydq, vi)
		u = H_var
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list