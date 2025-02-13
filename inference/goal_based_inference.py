#trying to implement the ideas in Lieberman & Willcox 2014

import sys
import math
import numpy as np
from mpmath import mp
import scipy.stats
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

sys.path.append('..')
	
#offline, train the model based on provided samples
def gbi_train_model(qoi_samples, y_samples, verbose=0, ncomp=0):
	#step 1: train 
	#Get a a Gaussian mixture model from the push-forward of the prior through yd and yp
	print("getting data...",flush=True) if verbose else 0
	yp_sample = qoi_samples
	yd_sample = y_samples

	p_dimension = 1 #len(yp_sample[0]) #1
	d_dimension = len(yd_sample[0]) #3

	#actually i think ill use sklearn for this
	data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_samples)])

	#need to use a "Bayesian information criterion approach [8]"
	#or, a "number of components to balance the maximization of likelihood of data and the Bayesian information criterion [4]"
	#which one??? the latter seems fine
	#ok, generally to do this, we'll search bottom-up, increasing n_comp until we find a minimal BIC
	gmm = None
	if ncomp == 0:
		print("calculating ideal ncomp...",flush=True) if verbose else 0
		curr_gmm = GaussianMixture(n_components=1).fit(data)
		curr_bic = curr_gmm.bic(data)
		print(curr_bic, flush=True) if verbose==2 else 0
		next_gmm = GaussianMixture(n_components=2).fit(data)
		next_bic = next_gmm.bic(data)
		print(next_bic, flush=True) if verbose==2 else 0
		ncomp = 1
		while next_bic < curr_bic:
			curr_bic = next_bic
			curr_gmm = next_gmm
			ncomp += 1
			next_gmm = GaussianMixture(n_components=ncomp+1).fit(data)
			next_bic = next_gmm.bic(data)
			print(next_bic, flush=True) if verbose else 0

		print("Number of components that minimizes BIC is",ncomp, flush=True) if verbose else 0
		gmm = curr_gmm
		
	else:
		#Allow someone to override this process if they think they know how many components they want
		gmm = GaussianMixture(n_components=ncomp).fit(data)
		print("Using provided number of components,",ncomp,flush=True) if verbose else 0
		
	return gmm

#online, condition the joint gmm on the data to get the posterior predictive
def gbi_condition_model(gmm, Yd, verbose=0):
	ncomp = gmm.n_components
	#step 2
	#Now we have our data, and we find the posterior predictive from that
	print("calculating posterior predictive...",flush=True) if verbose else 0
	
	#with mp.workdps(20):
	###
	Yd = mp.matrix(Yd) #column

	#get key parameters from the GMM
	mu = [mp.matrix(mu_k) for mu_k in gmm.means_]
	Sig = [mp.matrix(cov_k) for cov_k in gmm.covariances_]
	alpha = mp.matrix(gmm.weights_)
	p_dimension = len(mu[0]) - len(Yd) #scalar, this is usually 1
	ymean_p = [mu_k[:p_dimension] for mu_k in mu] #column
	ymean_d = [mu_k[p_dimension:] for mu_k in mu] #column (kinda)
	Sig_pp = [Sig_k[:p_dimension, :p_dimension] for Sig_k in Sig]
	Sig_pd = [Sig_k[:p_dimension, p_dimension:] for Sig_k in Sig]
	Sig_dp = [Sig_k[p_dimension:, :p_dimension] for Sig_k in Sig]
	Sig_dd = [Sig_k[p_dimension:, p_dimension:] for Sig_k in Sig]
	
	if verbose==2:
		print("ymean_p:\n", ymean_p)
		print("ymean_d:\n", ymean_d)
		print("Sig_pp:\n", Sig_pp)
		print("Sig_pd:\n", Sig_pd)
		print("Sig_dp:\n", Sig_dp)
		print("Sig_dd:\n", Sig_dd, flush=True)

	#parameters for the new GMM:
	B1 = [alpha[k] * (2*math.pi)**(-ncomp/2.0) * mp.det(Sig_dd[k])**(-0.5) 
			* mp.exp(mp.norm(-0.5 * (Yd - ymean_d[k]).T * (Sig_dd[k])**-1 * (Yd - ymean_d[k])),p=1)
			for k in range(ncomp)] #something is busted here
	B0 = sum(B1)
	if B0 == 0:
		print("Something stupid happened")
		print(B1)
		sys.exit()
	beta = [B1[k] / B0 for k in range(ncomp)]
	mu_Yd = [ymean_p[k] + Sig_pd[k] * (Sig_dd[k])**-1 * (Yd - ymean_d[k]) for k in range(ncomp)]
	Sig_Yd = [Sig_pp[k] - Sig_pd[k] * (Sig_dd[k])**-1 * Sig_dp[k] for k in range(ncomp)]
	###
	
	#convert back to np arrays, annoyingly:
	beta = np.array([float(mp.norm(x,p=1)) for x in beta])
	mu_Yd = np.array([[float(mp.norm(x,p=1))] for x in mu_Yd])
	Sig_Yd = np.array([[[float(mp.norm(x,p=1))]] for x in Sig_Yd])
	
	if verbose==2:
		print("beta:", beta)
		print("mu_Yd:", mu_Yd)
		print("Sig_Yd:", Sig_Yd, flush=True)
	
	return beta, mu_Yd, Sig_Yd
	
#online, get the pdf from the posterior predictive
#i.e., given the data Yd, what is the probability of seeing the QoI yp
def gbi_pdf_posterior_predictive(beta, mu_Yd, Sig_Yd, yp, verbose=0):
	#first, get the parameters of the posterior predictive
	beta, mu_Yd, Sig_Yd = gbi_condition_model(gmm, Yd, verbose=verbose)
	
	#I think I'm assuming here that yp is 1-dimensional, generalize this later...
	pdfs = np.array([p * scipy.stats.norm.pdf(x=yp, loc=mu, scale=np.sqrt(sd)) for mu, sd, p in zip(mu_Yd, Sig_Yd, beta)])
	pdf = np.sum(np.array(pdfs))

	return pdf
	
def gbi_gmm_variance(beta, mu_Yd, Sig_Yd):
	moment1 = np.sum([b*mu for b,mu in zip(beta, mu_Yd)])
	moment2 = np.sum([b*(var + mu**2) for b,mu,var in zip(beta, mu_Yd, Sig_Yd)])
	var = moment2 - moment1**2
	return var
	
def gbi_var_of_conditional_pp(gmm, Yd, verbose=0):
	beta, mu_Yd, Sig_Yd = gbi_condition_model(gmm, Yd, verbose=verbose)
	return gbi_gmm_variance(beta, mu_Yd, Sig_Yd)
	
#Given a GMM and a conditional Yd, what is a sample of Yp?
def gbi_sample_of_conditional_pp(gmm, Yd, verbose=0):
	beta, mu_Yd, Sig_Yd = gbi_condition_model(gmm, Yd, verbose=verbose)
	samples = np.array([p * scipy.stats.norm.rvs(loc=mu, scale=np.sqrt(sd)) for mu, sd, p in zip(mu_Yd, Sig_Yd, beta)])
	sample = np.sum(np.array(samples))
	
	return sample
	

def plot_predictive_posterior(beta, mu_Yd, Sig_Yd, lbound, rbound, drawplot=True, plotmean=False, compplot=True, maincolor='k'):
	#p is one-dimensional, give it a plot
	x = np.linspace(lbound, rbound, 10000)

	pdfs = np.array([p * scipy.stats.norm.pdf(x=x, loc=mu, scale=np.sqrt(sd)) for mu, sd, p in zip(mu_Yd, Sig_Yd, beta)])
	pdfs = np.array([pdf[0] for pdf in pdfs])
	density = np.sum(np.array(pdfs), axis=0)

	if compplot:
		for pdf in pdfs:
			plt.plot(x, pdf, '--', c='gray')#c=np.random.rand(3))
	plt.plot(x, density, '-', c=maincolor)
	plt.xlabel('$y_p$')
	plt.ylabel('$f(y_p | y_d)$')
	if plotmean:
		mean = np.sum([b*mu for b,mu in zip(beta,mu_Yd)])
		plt.axvline(mean, c='red')
	if drawplot:
		plt.show()
		
def find_ncomp(theta_samples, qoi_samples, y_samples):
	#step 1: train 
	#Get a a Gaussian mixture model from the push-forward of the prior through yd and yp
	samples = theta_samples
	yp_sample = qoi_samples
	yd_sample = y_samples

	p_dimension = 1 #len(yp_sample[0]) #1
	d_dimension = len(yd_sample[0]) #3

	#actually i think ill use sklearn for this
	data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(samples)])

	#need to use a "Bayesian information criterion approach [8]"
	curr_gmm = GaussianMixture(n_components=1).fit(data)
	curr_bic = curr_gmm.bic(data)
	next_gmm = GaussianMixture(n_components=2).fit(data)
	next_bic = next_gmm.bic(data)
	ncomp = 1
	while next_bic < curr_bic:
		curr_bic = next_bic
		curr_gmm = next_gmm
		ncomp += 1
		next_gmm = GaussianMixture(n_components=ncomp+1).fit(data)
		next_bic = next_gmm.bic(data)
		
	return ncomp


#train a model on the data p(theta),d -> y and p(theta) -> Q
#to develop the joint model Q = J(y,d)
#for use in the goal-based approach
#with the same GMM approach I used originally
def joint_model_gmm(problem, N, doPrint=False, doPlot=False):
	if doPrint:
			print("Generating the training data...", flush=True)
			
	#model inputs
	theta_train = problem.prior_rvs(N)
	d_train = [d for d in problem.sample_d(N)]
	
	#model outputs
	qoi_train = [problem.H(theta) for theta in theta_train]
	y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
	
	#print(theta_train)
	#print(d_train)
	#print(qoi_train)
	#print(y_train)
	
	joint_input = [yi+di for yi,di in zip(y_train, d_train)]
	joint_output = np.array(qoi_train)
	
	#print(joint_input)
	#print(joint_output)
	
	if doPrint:
			print("Training the joint model...", flush=True)

	# Create linear regression object
	gmm = gbi_train_model(joint_output, joint_output, joint_input, ncomp=0, verbose=doPrint)
	
	model = gmm
	
	return model