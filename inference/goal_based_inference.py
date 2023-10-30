#trying to implement the ideas in Lieberman & Willcos 2014

import sys
import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

sys.path.append('..')
	
#offline, train the model based on provided samples
def gbi_train_model(theta_samples, qoi_samples, y_samples, verbose=0):
	#step 1: train 
	#Get a a Gaussian mixture model from the push-forward of the prior through yd and yp
	print("getting data...",flush=True) if verbose else 0
	samples = theta_samples
	yp_sample = qoi_samples
	yd_sample = y_samples

	p_dimension = 1 #len(yp_sample[0]) #1
	d_dimension = len(yd_sample[0]) #3

	#actually i think ill use sklearn for this
	data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(samples)])

	#need to use a "Bayesian information criterion approach [8]"
	#or, a "number of components to balance the maximization of likelihood of data and the Bayesian information criterion [4]"
	#which one??? the latter seems fine
	#ok, generally to do this, we'll search bottom-up, increasing n_comp until we find a minimal BIC
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
	return gmm

#online, condition the joint gmm on the data to get the posterior predictive
def gbi_condition_model(gmm, Yd, verbose=0):
	ncomp = gmm.n_components
	#step 2
	#Now we have our data, and we find the posterior predictive from that
	print("calculating posterior predictive...",flush=True) if verbose else 0
	Yd = np.transpose(np.array(Yd))

	#get key parameters from the GMM
	mu = np.array(gmm.means_)
	Sig = np.array(gmm.covariances_)
	alpha = np.array(gmm.weights_)
	p_dimension = len(mu[0]) - len(Yd) #this is usually 1
	ymean_p = np.array([mu_k[:p_dimension] for mu_k in mu])
	ymean_d = np.array([np.transpose(mu_k[p_dimension:]) for mu_k in mu])
	Sig_pp = np.array([Sig_k[:p_dimension, :p_dimension] for Sig_k in Sig])
	Sig_pd = np.array([Sig_k[:p_dimension, p_dimension:] for Sig_k in Sig])
	Sig_dp = np.array([Sig_k[p_dimension:, :p_dimension] for Sig_k in Sig])
	Sig_dd = np.array([Sig_k[p_dimension:, p_dimension:] for Sig_k in Sig])
	
	if verbose==2:
		print("ymean_p:\n", ymean_p)
		print("ymean_d:\n", ymean_d)
		print("Sig_pp:\n", Sig_pp)
		print("Sig_pd:\n", Sig_pd)
		print("Sig_dp:\n", Sig_dp)
		print("Sig_dd:\n", Sig_dd, flush=True)

	#parameters for the new GMM:
	B1 = [alpha[k] * (2*math.pi)**(-ncomp/2.0) * np.linalg.det(Sig_dd[k])**(-0.5) 
			* np.exp(-0.5 * np.transpose(Yd - ymean_d[k]) @ np.linalg.inv(Sig_dd[k]) @ (Yd - ymean_d[k])) 
			for k in range(ncomp)]
	B0 = sum(B1)
	beta = np.array([B1[k] / B0 for k in range(ncomp)])
	mu_Yd = np.array([ymean_p[k] + Sig_pd[k] @ np.linalg.inv(Sig_dd[k]) @ (Yd - ymean_d[k]) for k in range(ncomp)])
	Sig_Yd = np.array([Sig_pp[k] - Sig_pd[k] @ np.linalg.inv(Sig_dd[k]) @ Sig_dp[k] for k in range(ncomp)])
	
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

def plot_predictive_posterior(beta, mu_Yd, Sig_Yd, lbound, rbound, drawplot=True):
	#p is one-dimensional, give it a plot
	x = np.linspace(lbound, rbound, 10000)

	pdfs = np.array([p * scipy.stats.norm.pdf(x=x, loc=mu, scale=np.sqrt(sd)) for mu, sd, p in zip(mu_Yd, Sig_Yd, beta)])
	pdfs = np.array([pdf[0] for pdf in pdfs])
	density = np.sum(np.array(pdfs), axis=0)

	for pdf in pdfs:
		plt.plot(x, pdf, '--', c='gray')#c=np.random.rand(3))
	plt.plot(x, density, '-k')
	plt.xlabel('$y_p$')
	plt.ylabel('$f(y_p | y_d)$')
	if drawplot:
		plt.show()