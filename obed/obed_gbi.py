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
from inference.bn_modeling import *

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
	
def U_varH_gbi_joint(d, problem, gmm_qyd, n_mc, doPrint=False):   
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	if doPrint:
		print("Generating",n_mc,"joint MC samples of theta and y...",flush=True)
	pthetas = problem.prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	if doPrint:
		print(Y1_list, flush=True)
	
	#We expect we have already trained a joint gmm model p(y,d,q) offline,
	#and will condition it for p(q|y,d)
	if doPrint:
		print("Conditioning GMM in MC loop...",flush=True)
		
	#Pre-calculate an inverse matrix, to speed up the MC loop:
	inv_Sig_dd, logdet_Sig_dd = gbi_precalc_Sigdd(gmm_qyd, p_dim=1)
	
	U_list = []
	for i,y in enumerate(Y1_list): #MC loop		
		vi = np.concatenate((y,d))
	
		#Now, use my posterior predictive to calculate the utility
		H_var = gbi_var_of_conditional_pp(gmm_qyd, vi, 
			inv_Sig_dd_precalc=inv_Sig_dd, logdet_Sig_dd_precalc=logdet_Sig_dd)
		u = H_var
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list

#Like U_varH_gbi_joint, except instead of drawing new y and theta samples,
#we draw a randomized n_mc subset from a provided qyd sample file,
#and feed it in here
def U_varH_gbi_joint_presampled(d, problem, gmm_qyd, presampled_ylist, n_mc, doPrint=False):   
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	if doPrint:
		print("Pulling",n_mc,"presampled joint MC samples of theta and y...",flush=True)
	mc_ylist = bn_random_subset(presampled_ylist, n_mc, allowDuplicates=True)
	
	#We expect we have already trained a joint gmm model p(y,d,q) offline,
	#and will condition it for p(q|y,d)
	if doPrint:
		print("Conditioning GMM in MC loop...",flush=True)
		
	#Pre-calculate an inverse matrix, to speed up the MC loop:
	inv_Sig_dd, logdet_Sig_dd = gbi_precalc_Sigdd(gmm_qyd, p_dim=1)
	
	U_list = []
	for i,y in enumerate(mc_ylist): #MC loop		
		vi = np.concatenate((y,d))
	
		#Now, use my posterior predictive to calculate the utility
		H_var = gbi_var_of_conditional_pp(gmm_qyd, vi, inv_Sig_dd_precalc=inv_Sig_dd, 
			logdet_Sig_dd_precalc=logdet_Sig_dd)
		u = H_var
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list
	
#information criterion
def U_infoH_gbi_joint(d, problem, gmm_qyd, gmm_q, n_mc, doPrint=False):   
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	if doPrint:
		print("Pulling",n_mc,"samples of theta_i, propagated to y_i and Q_i...",flush=True)
	theta_ilist = problem.prior_rvs(n_mc)
	y_ilist = [problem.eta(theta, d) for theta in theta_ilist]
	Q_ilist = [problem.H(theta) for theta in theta_ilist]
	
	#We expect we have already trained a joint gmm model p(y,d,q) offline,
	#and will condition it for p(q|y,d)
	if doPrint:
		print("Conditioning GMM in MC loop...",flush=True)
		
	#Pre-calculate an inverse matrix, to speed up the MC loop:
	inv_Sig_dd, logdet_Sig_dd = gbi_precalc_Sigdd(gmm_qyd, p_dim=1)
	
	#Grab parameters from the GMM of Q
	beta_Q = gmm_q.weights_
	mu_Q = gmm_q.means_
	Sig_Q = gmm_q.covariances_
	print(len(beta_Q), len(mu_Q), len(Sig_Q))
	#TODO Filter out all of the zero weight here
	
	#plot_predictive_posterior(beta_Q, mu_Q, Sig_Q, 0, 5, drawplot=True, plotmean=False, compplot=True, maincolor='k')
	
	U_list = []
	for i,yi in enumerate(y_ilist): #MC loop		
		vi = np.concatenate((yi,d))
		Qi = Q_ilist[i]
	
		#First, condition the joint GMM to get a GMM on Q|y,d
		beta_cond, mu_cond, Sig_cond = gbi_condition_model(gmm_qyd, vi, inv_Sig_dd_precalc=inv_Sig_dd, logdet_Sig_dd_precalc=logdet_Sig_dd, verbose=False)		
		#plot_predictive_posterior(beta_cond, mu_cond, Sig_cond, -100, 100, drawplot=True, plotmean=False, compplot=True, maincolor='k')
		
		print(len(beta_cond), len(mu_cond), len(Sig_cond))
		#De-normalize the GMMs here? Or can i find a way to do it outside the loop?
		#Actually no, gbi_condition_model returns de-normalized GMM parameters
		#And I can just train gmm_q to never have been normalized
		
		#TODO Filter out all of the zero weight here
		
		#Then, solve DKL[ p(Q|y,d) || p(Q) ]
		X = [Qi] #single sample, single feature
		#TODO this is broken for zero weight. Need to filter them and do it manually
		#q_yd_X = gmm_q_yd.score_samples(X) #computes log likelihood of X
		#q_X = gmm_q.score_samples(X) #computes log likelihood of X
		#add all weighted scores, then take the log of that sum
		q_yd_X = sum([w*scipy.stats.norm.logpdf(X, loc=mu_cond[i], scale=np.sqrt(Sig_cond[i])) for i,w in enumerate(beta_cond)]).item()
		q_X = sum([w*scipy.stats.norm.logpdf(X, loc=mu_Q[i], scale=np.sqrt(Sig_Q[i])) for i,w in enumerate(beta_Q)]).item()
		#print(q_yd_X,q_X)
		u = np.log(q_yd_X) - np.log(q_X)
		
		#TODO handle log errors
			
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list