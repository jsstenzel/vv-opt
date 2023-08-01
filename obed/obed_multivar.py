#This details an obed script that can be called regardless of the problem
#so it takes as input some function defintions

#MCMC References:
#https://twiecki.io/blog/2015/11/10/mcmc-sampling/
#https://towardsdatascience.com/bayesian-inference-and-markov-chain-monte-carlo-sampling-in-python-bada1beabca7

#Not a code reference, but provides a nice example for how to tune and fiddle mcmc params:
#https://github.com/gmcgoldr/pymcmc

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt
import seaborn as sns


		

def multivar_likelihood_kernel(d, exp_fn, prior_rvs, n1):
	thetas = prior_rvs(n1)
	Y1_list = [exp_fn(theta, d) for theta in thetas]
	
	y_theta_values = np.array([np.concatenate([thetas[i],Y1_list[i]]) for i,_ in enumerate(thetas)])
	likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)

	return likelihood_kernel, Y1_list
	
	
#if you want to set the ythetas manually, use this
def general_likelihood_kernel(*params):
	#first ensure all lists have equal length
	if False in [len(i) == len(params[0]) for i in params]:
		print("general_likelihood_kernel failure: all inputs must have same length,", [len(i) for i in params])
		
	#then put lists together
	y_theta_values = np.array([np.concatenate([param[i] for param in params]) for i,_ in enumerate(params[0])])

	likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)			
	return likelihood_kernel

#likelihood is the kde generated from sample
#this assumes 3 elements in 7 - i can generalize this easily later
def likelihood_plot(likelihood, sample):
	xdata = [y[0] for y in sample]
	ydata = [y[1] for y in sample]
	zdata = [y[2] for y in sample]
	
	#3d plot
	#fig = plt.figure()
	#ax = plt.axes(projection='3d')
	#ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')
	
	xs = np.linspace(min(xdata), max(xdata), len(xdata))
	plt.plot(xs, [likelihood([x,np.mean(ydata),np.mean(zdata)]) for x in xs]    )
	plt.show()
	
	ys = np.linspace(min(ydata), max(ydata), len(ydata))
	plt.plot(ys, [likelihood([np.mean(xdata),y,np.mean(zdata)]) for y in ys]    )
	plt.show()
	
	zs = np.linspace(min(zdata), max(zdata), len(zdata))
	plt.plot(zs, [likelihood([np.mean(xdata),np.mean(ydata),z]) for z in zs]    )
	plt.show()

#y                 - single data point
#likelihood_kernel - function that you can evaluate to determine the likelihood fn at theta, y
#proposal_dist     - distribution fn used to propose new data points; takes a mean as argument
def mcmc_multivar(y, likelihood_kernel, proposal_dist, prior_rvs, prior_pdf_unnorm, n, burnin=0, lag=1, doPlot=False, legend=None):
	N = n*lag + burnin #number of samples of the posterior i want, times lag plus burn-in
	#I will probably have to take into account things like burn-in/convergence and lag? idk
	
	theta_current = prior_rvs(1)
	mcmc_trace = []
	for i in range(N): #MCMC loop
		#Propose a new value of theta with Markov chain
		#The key here is that we want to generate a new theta randomly, and only dependent on the previous theta
		#Usual approach is a gaussian centered on theta_current with some efficient proposal_width, but that doesnt work on all domains
		#I'll keep it general: have an proposal distribution as an argument, which takes theta_current as a mean/argument itself
		theta_proposal = proposal_dist(theta_current)
		
		#Compute likelihood of the "data" for both thetas - reusing the loop-1 likelihood evaluations
		ytheta_current = np.concatenate([theta_current, y])
		ytheta_proposal = np.concatenate([theta_proposal, y])
		loglikelihood_current = np.log(likelihood_kernel(ytheta_current.T)[0]) #using kde to estimate pdf
		loglikelihood_proposal = np.log(likelihood_kernel(ytheta_proposal.T)[0]) #calc probability at the data y
		#Compute acceptance ratio
		logprior_current = np.log(np.prod(prior_pdf_unnorm(theta_current)))
		logprior_proposal = np.log(np.prod(prior_pdf_unnorm(theta_proposal)))
		logp_current = loglikelihood_current + logprior_current
		logp_proposal = loglikelihood_proposal + logprior_proposal
		logR = logp_proposal - logp_current
		R = np.exp(logR)
		#Accept our new value of theta according to the acceptance probability R:
		if np.random.random_sample() < R:
			theta_current = theta_proposal
		#Include theta_current in my trace according to rules
		if i > burnin and i%lag == 0:
			mcmc_trace.append(theta_current)
			
	if doPlot:
		plt.plot(mcmc_trace)
		plt.legend(legend)
		plt.show()
	return mcmc_trace


"""
This function is part of the OBED problem, solving the utility for a given d
This may be an inefficient brute-force approach
This assumes a utility u(d,y,theta) = 1/Var[H(theta|y,d)], minimizing the variance of the HLVA evaluated over the posterior
A higher function must optimize over the outputs of this model to find d*, the d that maximized U_d
"""
#This is like U_brute_varH, except we're reusing loop 1 samples in loop 3, like Huan & Marzouk
#Here, we're using straightforward Metropolis-Hastings to solve the multivariate MCMC problem
#This calculates the probability of meeting the requirement in distribution - used as constraint for the cost-optimization
def U_probreq(d, problem, mcmc_proposal, theta_domains, minreq=None, maxreq=None, n_mc=10**4, n_mcmc=10**3, n_kde=10**4, burnin=0, lag=1):   
	#requirement: do min, max, or both, but not neither
	if minreq==None and maxreq==None:
		print("U_probreq needs a requirement on the QoI in order to determine probability")
		exit()
	
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	pthetas = prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	
	#Then i want to create an estimated pdf of y as a fn of theta
	#I'll generate thetas uniformly across the domain (theta_domains), and sample ys from that
	kde_thetas = [[scipy.stats.uniform.rvs(size=1, loc=left, scale=right-left)[0] for left,right in theta_domains] for _ in range(n_kde)]
	kde_ys = [fp.eta(theta, d_historical) for theta in kde_thetas]
	likelihood_kernel = general_likelihood_kernel(kde_thetas, kde_ys)
	
	U_list = []
	for y in Y1_list: #MC loop		
		mcmc_trace = mcmc_multivar(y, likelihood_kernel, mcmc_proposal, problem.prior_rvs, problem.prior_pdf_unnorm, n_mcmc, burnin, lag)
			
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
		
	#compute an in-distribution probability
	U = np.average(U_list)
	#is mean the best thing for this? For most problems we'll be bounded by the high end much more than the low...
	#itll be interesting to see the U_list i guess
	return U, U_list

