#This details an obed script that can be called regardless of the problem

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
import random

sys.path.append('..')
from obed.mcmc import *

"""
this is a general evaluation of the utiltiy
assumes nothing about how we sample from the likelihood fn, or evaluate the likelihood of data
"""
def U_Huan(d, problem, G, n_out, n_in=0, doReuse=False, doPrint=False):   
	#Draw theta_i from the prior
	thetas = problem.prior_rvs(n_out)
	#these can be reused, potentially
	G_evals = [G(theta,d) for theta in thetas]
	#and draw y_i from theta_i and the likelihood
	ys = [problem.eta(theta_i, d, G_evals[i]) for i,theta_i in enumerate(thetas)]
		
	likelihood_list = []
	evidence_list = []
	U_list = []
	skipcount = 0
	for i,yi in enumerate(ys): #MC loop		
		#get the inner thetas for the evidence loop
		if n_in == 0:
			theta_js = thetas[:]
			G_js = G_evals[:]
		elif n_in>0 and doReuse==False:
			theta_js = problem.prior_rvs(n_in)
			G_js = [G(theta_j,d) for theta_j in theta_js]
		elif doReuse==True and n_in <= n_out:
			theta_js = thetas[:n_in]
			G_js = G_evals[:n_in]
		elif doReuse==True and n_in > n_out:
			#reuse what you can, and resample the rest
			thetaj_extra = problem.prior_rvs(n_in-n_out)
			theta_js = np.concatenate([thetas[:], thetaj_extra])
			G_js = np.concatenate([G_evals[:], [G(theta_j,d) for theta_j in thetaj_extra ]] )
		elif doReuse==2 and n_in <= n_out:
			rand_indices = [r for r in range(n_out)]
			random.shuffle(rand_indices)
			rand_indices = rand_indices[:n_in]
			theta_js = [thetas[r] for r in rand_indices]
			G_js = [G_evals[r] for r in rand_indices]
		else:
			print("error")
			quit()
		
		#Calculate the evidence, a list of p(yi|d) for yi in ys
		evidences_j = [problem.eta_likelihood(yi, theta_j, d, G_js[j]) for j,theta_j in enumerate(theta_js)]
		evidence_yi = np.mean(evidences_j)
		
		#evaluate the likelihood
		theta_i = thetas[i]
		likelihood = problem.eta_likelihood(yi, theta_i, d, G_evals[i])
	
		#guard against near-zero evidence estimates and infinite utilities:
		if evidence_yi < 1e-12:
			skipcount+=1
			continue
	
		#Now, calculate the utility
		u = np.log(likelihood) - np.log(evidence_yi)
		U_list.append(u)
		likelihood_list.append(likelihood)
		evidence_list.append(evidence_yi)
		if doPrint:
			print(str(i+1)+"/"+str(n_out),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.mean(U_list)
	return U, U_list, likelihood_list, evidence_list, skipcount

"""
Using no approximation for eta
We solve the likelihood fn by adding a 3rd Monte Carlo loop, and making a Gaussian kde over y for fixed theta
"""
def U_Huan_direct_3loop(d, problem, n_out, n_3rd, n_in=0, doPrint=False):   
	#Draw theta_i from the prior, and y_i from the likelihood
	thetas = problem.prior_rvs(n_out)
	ys = [problem.eta(theta_i, d) for theta_i in thetas]
		
	U_list = []
	for i,yi in enumerate(ys): #MC loop		
		#get the inner thetas for the evidence loop
		if n_in>0:
			theta_js = problem.prior_rvs(n_in)
			y_js = [problem.eta(theta_j, d) for theta_j in theta_js]
		else: #do reuse!
			theta_js = thetas[:]
			y_js = ys[:]
		
		#Calculate the evidence, a list of p(yi|d) for yi in ys
		evidence_yi = np.mean(y_js)
		
		#evaluate the likelihood
		theta_i = thetas[i]
		#first, generate eta(thetai,d) samples
		y_likelihood_samples = np.array([problem.eta(theta_i, d) for _ in range(n_3rd)])
		#then, generate a GMM with this sample
		likelihood_kde = scipy.stats.gaussian_kde(y_likelihood_samples.T)
		#lastly, evaluate the likelihood of the data yi
		likelihood = likelihood_kde(yi)
	
		#Now, calculate the utility
		u = np.log(likelihood) - np.log(evidence_yi)
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_out),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.mean(U_list)
	return U, U_list
	
"""
Using no approximation for eta
We solve the likelihood fn by adding another loop at level 2, to develop a Gaussian kde over theta and y
Can reuse for this loop as well
"""
def U_Huan_direct_parallelloop(d, problem, n_out, n_in=0, n_3rd=0, doPrint=False):   
	#Draw theta_i from the prior, and y_i from the likelihood
	thetas = problem.prior_rvs(n_out)
	ys = [problem.eta(theta_i, d) for theta_i in thetas]
	
	###outside the main loop, develop a kde for y and theta
	if n_3rd>0:
		theta_likelihood = problem.prior_rvs(n_3rd)
		y_likelihood = [problem.eta(theta, d) for theta in theta_likelihood]
	else: #do reuse!
		theta_likelihood = thetas[:]
		y_likelihood = ys[:]
	#first, generate eta(thetai,d) samples
	y_theta_values = np.array([np.concatenate([theta_likelihood[i],y_likelihood[i]]) for i,_ in enumerate(theta_likelihood)])
	#then, generate a GMM with this sample
	likelihood_kde = scipy.stats.gaussian_kde(y_theta_values.T)
	
	U_list = []
	for i,yi in enumerate(ys): #MC loop		
		#get the inner thetas for the evidence loop
		if n_in>0:
			theta_js = problem.prior_rvs(n_in)
			y_js = [problem.eta(theta_j, d) for theta_j in theta_js]
		else: #do reuse!
			theta_js = thetas[:]
			y_js = ys[:]
		
		#Calculate the evidence, a list of p(yi|d) for yi in ys
		evidence_yi = np.mean(y_js)
		
		#evaluate the likelihood of the data yi
		theta_i = thetas[i]
		ytheta_i = np.concatenate([theta_i,yi])
		likelihood = likelihood_kde(ytheta_i)
	
		#Now, calculate the utility
		u = np.log(likelihood) - np.log(evidence_yi)
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_out),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.mean(U_list)
	return U, U_list




def U_Huan_PCE(d, problem, eta_polynomial, n_out, n_in=0, doReuse=True, doPrint=False): 
	0