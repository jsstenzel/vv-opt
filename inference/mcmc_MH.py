#MCMC References:
#https://twiecki.io/blog/2015/11/10/mcmc-sampling/
#https://towardsdatascience.com/bayesian-inference-and-markov-chain-monte-carlo-sampling-in-python-bada1beabca7

#Not a code reference, but provides a nice example for how to tune and fiddle mcmc params:
#https://github.com/gmcgoldr/pymcmc

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acf


sys.path.append('..')
from obed.pdf_estimation import *

def prop_fn(theta_prev, prop_width):
	return scipy.stats.multivariate_normal.rvs(mean=theta_prev, cov=prop_width, size=1)

#TODO i need to hack this to make it very general
def mcmc_MH(n_mcmc, n_trunc, prop_Sigma, doAdaptive=False, doPrint=False):	
	#get data
	f = open('inferencedata.json')
	data = json.load(f)
	x_pts = data['xobserved']
	y_data = data['Uobserved']
	
	#TODO
	# Optimal scaling factor for the proposal covariance (Gelman et al.)
	sd = (2.38 ** 2) / d
	
	Z_current = scipy.stats.norm.rvs(size=n_trunc) #sample from prior
	likelihoodfn_current = full_stochastic_likelihood(x_pts, y_data, Z_current)
	rolling_Zmean = Z_current
	trace_covariance = prop_Sigma
	
	def prop_fn(Z_prev, prop_width):
		return scipy.stats.multivariate_normal.rvs(mean=Z_prev, cov=prop_width, size=1)
	
	mcmc_trace = []
	acceptance_count = 0
	randwalk_count = 0
	for i in range(n_mcmc): #MCMC loop
		#Propose a new value of theta with Markov chain
		Z_proposal = prop_fn(Z_current, prop_Sigma)
		
		#Compute likelihood of the data y for both thetas
		likelihoodfn_proposal = full_stochastic_likelihood(x_pts, y_data, Z_proposal)

		#Compute acceptance ratio
		prior_current = np.prod(scipy.stats.norm.pdf(Z_current))
		prior_proposal = np.prod(scipy.stats.norm.pdf(Z_proposal))
		
		R = (prior_proposal * likelihoodfn_proposal) / (prior_current * likelihoodfn_current)
		#print(R)

		#Accept our new value of theta according to the acceptance probability R:
		if np.random.random_sample() < R:
			#accept
			Z_current = Z_proposal
			likelihoodfn_current = likelihoodfn_proposal
			acceptance_count += 1
			
		acceptance_rate = acceptance_count / (i+1)
		mcmc_trace.append(Z_current.tolist())
		
		"""
		if doAdaptive:
			ii=i+1
			if ii > 1:
				rolling_Zmean = (ii/(ii+1)) * rolling_Zmean + (1/(ii+1)) * Z_current
				delta_n = np.array(Z_current) - np.array(rolling_Zmean)
				trace_covariance = trace_covariance * ((ii-2)/(ii-1)) + (1/ii) * (delta_n @ delta_n.T)
			if np.random.random_sample() < R:
					prop_Sigma = trace_covariance
		"""
		if doAdaptive and (i>200):
			if np.random.random_sample() < R: #only update the covariance upon acceptance?
				data = np.transpose(mcmc_trace)
				covariance = np.cov(data, ddof=1)
				prop_Sigma = covariance
			
		if doPrint:
			print(str(i+1)+"/"+str(n_mcmc), acceptance_count, '\t',  flush=True, end='\r')
			

	acceptance_rate = acceptance_count / n_mcmc
	return mcmc_trace, acceptance_rate

	
def autocorrelate_mcmc(trace, burnin=0, lag=1, doPlot=True):
	tr = np.array(trace[burnin::lag])
	data = np.transpose(tr)
	
	for i,zi in enumerate(data):
		#autocorr = np.correlate(zi, zi, mode='same')
		autocorr = acf(zi, nlags=len(tr)-2)
		#print(autocorr)
		plt.plot(autocorr)
		plt.title("MCMC autocorrelation for Z_"+str(i))
		plt.xlabel("i")
		plt.show()

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
			#plt.xscale('log')
			plt.plot(range(1,len(data_i)+1), data_i, c='r')
			#plt.plot([np.mean(data_i[0:n]) for n in range(len(data_i))], c='g')
			#plt.plot([np.std(data_i[0:n], ddof=1) for n in range(len(data_i))], c='b')
			#plt.legend(["trace","mean trace","stddev trace"])
			plt.xlabel("n runs")
			plt.title("MCMC Trace for theta["+str(i)+"]")
			plt.show()
	
	return means, stddevs, covariance
