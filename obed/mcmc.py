#MCMC References:
#https://twiecki.io/blog/2015/11/10/mcmc-sampling/
#https://towardsdatascience.com/bayesian-inference-and-markov-chain-monte-carlo-sampling-in-python-bada1beabca7

#Not a code reference, but provides a nice example for how to tune and fiddle mcmc params:
#https://github.com/gmcgoldr/pymcmc

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

sys.path.append('..')
from obed.pdf_estimation import *

#y				 - single data point
#likelihood_kernel - function that you can evaluate to determine the likelihood fn at theta, y
#prop_fn	 - distribution fn used to propose new data points; takes a mean as argument
def mcmc_kernel(y, likelihood_kernel, prop_fn, prop_width, prior_rvs, prior_pdf_unnorm, n, burnin=0, lag=1, doPlot=False, legend=None):
	N = n*lag + burnin #number of samples of the posterior i want, times lag plus burn-in
	
	theta_current = prior_rvs(1)
	mcmc_trace = []
	acceptance_count = 0
	randwalk_count = 0
	for i in range(N): #MCMC loop
		#Propose a new value of theta with Markov chain
		#The key here is that we want to generate a new theta randomly, and only dependent on the previous theta
		#Usual approach is a gaussian centered on theta_current with some efficient proposal_width, but that doesnt work on all domains
		#I'll keep it general: have an proposal distribution as an argument, which takes theta_current as a mean/argument itself
		theta_proposal = prop_fn(theta_current, prop_width)
		
		#Compute likelihood of the "data" for both thetas with the provided kernel
		ytheta_current = np.concatenate([theta_current, y])
		ytheta_proposal = np.concatenate([theta_proposal, y])
		loglikelihood_current = likelihood_kernel.logpdf(ytheta_current.T)[0] #using kde to estimate pdf
		loglikelihood_proposal = likelihood_kernel.logpdf(ytheta_proposal.T)[0] #calc probability at the data y
		
		#Compute acceptance ratio
		prior_current = np.prod(prior_pdf_unnorm(theta_current))
		logp_current = loglikelihood_current + np.log(prior_current)
		if(not np.isfinite(logp_current) or prior_current==0):
			#Crazy idea: for R=nan, random walk until you get to suitable region
			R = 1 #dont even check p_proposal, just accept it
			randwalk_count += 1
		else:
			prior_proposal = np.prod(prior_pdf_unnorm(theta_proposal))
			if prior_proposal==0:
				R=0
			else:
				logp_proposal = loglikelihood_proposal + np.log(prior_proposal)
				logR = logp_proposal - logp_current
				R = np.exp(logR)
		#print(R)
		
		#Accept our new value of theta according to the acceptance probability R:
		if np.random.random_sample() < R:
			theta_current = theta_proposal
			acceptance_count += 1
		#Include theta_current in my trace according to rules
		if i > burnin and i%lag == 0:
			mcmc_trace.append(theta_current)
			
		acceptance_rate = acceptance_count / (i+1)
		randwalk_rate = randwalk_count / (i+1)
		print(str(i+1)+"/"+str(N), acceptance_count, randwalk_count, theta_current, '\t', end='\r', flush=True)
		
	acceptance_rate = acceptance_count / N
	randwalk_rate = randwalk_count / N
		
	if doPlot:
		plt.plot(mcmc_trace)
		plt.legend(legend)
		plt.show()
	return mcmc_trace, acceptance_rate, randwalk_rate


#y				 - single data point
#proposal_dist	 - distribution fn used to propose new data points; takes a mean as argument
#eta			   - model function which samples from the conditional distribution p(y|theta,d)
def mcmc_nolikelihood(y, d, eta, prop_fn, prop_width, prior_rvs, prior_pdf_unnorm, n_mcmc, n_kde, burnin=0, lag=1, doPlot=False, legend=None):
	N = n_mcmc*lag + burnin #number of samples of the posterior i want, times lag plus burn-in
	
	theta_current = prior_rvs(1)
	#need to estimate probability dists p(y|theta=theta_current,d=d)
	ysample_current = [eta(theta_current, d) for _ in range(n_kde)]
	#print("ycurr",np.mean(ysample_current, axis=0))
	likelihoodfn_current = general_likelihood_kernel(ysample_current)
	
	mcmc_trace = []
	acceptance_count = 0
	randwalk_count = 0
	for i in range(N): #MCMC loop
		#Propose a new value of theta with Markov chain
		#The key here is that we want to generate a new theta randomly, and only dependent on the previous theta
		#Usual approach is a gaussian centered on theta_current with some efficient proposal_width, but that doesnt work on all domains
		#I'll keep it general: have an proposal distribution as an argument, which takes theta_current as a mean/argument itself
		theta_proposal = prop_fn(theta_current, prop_width)
		
		#Compute likelihood of the "data" y for both thetas
		#So i need to estimate probability dist p(y|theta=theta_proposal,d=d)
		ysample_proposal = [eta(theta_proposal, d) for _ in range(n_kde)]
		#print("yprop",np.mean(ysample_proposal, axis=0))
		likelihoodfn_proposal = general_likelihood_kernel(ysample_proposal)

		loglikelihood_current = likelihoodfn_current.logpdf(y)[0]
		loglikelihood_proposal = likelihoodfn_proposal.logpdf(y)[0]
		#Compute acceptance ratio
		prior_current = np.prod(prior_pdf_unnorm(theta_current))
		logp_current = loglikelihood_current + np.log(prior_current)
		if(not np.isfinite(logp_current) or prior_current==0):
			#Crazy idea: for R=nan, random walk until you get to suitable region
			R = 1 #dont even check p_proposal, just accept it
			randwalk_count += 1
		else:
			prior_proposal = np.prod(prior_pdf_unnorm(theta_proposal))
			if prior_proposal==0:
				R=0
			else:
				logp_proposal = loglikelihood_proposal + np.log(prior_proposal)
				logR = logp_proposal - logp_current
				R = np.exp(logR)
		#print(R)

		#Accept our new value of theta according to the acceptance probability R:
		if np.random.random_sample() < R:
			theta_current = theta_proposal
			likeihoodfn_current = likelihoodfn_proposal
			acceptance_count += 1
		#Include theta_current in my trace according to rules
		if i > burnin and i%lag == 0:
			mcmc_trace.append(theta_current)
			
		acceptance_rate = acceptance_count / (i+1)
		randwalk_rate = randwalk_count / (i+1)
		print(str(i+1)+"/"+str(N), acceptance_count, randwalk_count, theta_current, '\t', end='\r', flush=True)
			
	acceptance_rate = acceptance_count / N
	randwalk_rate = randwalk_count / N
	
	if doPlot:
		plt.plot(mcmc_trace)
		plt.legend(legend)
		plt.show()
	return mcmc_trace, acceptance_rate, randwalk_rate


def mcmc_analyze(trace, burnin=0, lag=1, doPlot=False):
	#optionally perform more burnin and lag here
	#a[start_index:end_index:step]
	tr = np.array(trace[burnin::lag])
	data = np.transpose(tr)
	covariance = np.cov(data, ddof=1)
	
	means = []
	stddevs = []	
	
	for i,data_i in enumerate(data):
		means.append(np.mean(data_i))
		stddevs.append(np.std(data_i, ddof=1))
		if doPlot: #plotting the trace of the data, the mean, and the stddev
			plt.xscale('log')
			plt.plot(range(1,len(data_i)+1), data_i, c='r')
			plt.plot([np.mean(data_i[0:n]) for n in range(len(data_i))], c='g')
			plt.plot([np.std(data_i[0:n], ddof=1) for n in range(len(data_i))], c='b')
			plt.legend(["trace","mean trace","stddev trace"])
			plt.xlabel("n runs")
			plt.title("MCMC Trace for theta["+str(i)+"]")
			plt.show()
	
	return means, stddevs, covariance