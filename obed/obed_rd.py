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
	
def U_rd_info_gbi(d, problem, gmm_qyd, n_mc, doPrint=False):   
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
		#H_var = gbi_var_of_conditional_pp(gmm_qyd, vi, 
		#	inv_Sig_dd_precalc=inv_Sig_dd, logdet_Sig_dd_precalc=logdet_Sig_dd)
		#u = H_var
		u=0
		
		#TODO use conditioned gmm to find posterior of q
		#TODO calculate the criterion
		
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list