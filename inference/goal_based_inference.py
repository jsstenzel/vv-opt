#trying to implement the ideas in Lieberman & Willcox 2014

import sys
import math
import numpy as np
from mpmath import mp
import scipy.stats
from scipy.special import logsumexp
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

sys.path.append('..')

class GaussianMixtureNormalized(GaussianMixture):
	def __init__(self, n_components=1, *, covariance_type='full', tol=0.001, reg_covar=1e-06, max_iter=100, n_init=1, init_params='kmeans', weights_init=None, means_init=None, precisions_init=None, random_state=None, warm_start=False, verbose=0, verbose_interval=10, standardized_mean, standardized_std):
		self.standardized_mean = standardized_mean
		self.standardized_std = standardized_std
		super().__init__(n_components=n_components, covariance_type=covariance_type, tol=tol, reg_covar=reg_covar, max_iter=max_iter, n_init=n_init, init_params=init_params, weights_init=weights_init, means_init=means_init, precisions_init=precisions_init, random_state=random_state, warm_start=warm_start, verbose=verbose, verbose_interval=verbose_interval)
	
	"""	
	#I dont use this anywhere, but it would be necessary to do this for completeness
	def sample(n_samples=1):
		samples = super().sample(n_samples)
		destandardized_samples = [] #TODO
		return destandardized_samples
	"""

#offline, train the model based on provided samples
def gbi_train_model(qoi_samples, y_samples, verbose=0, ncomp=0, ncomp_start=1, careful=False):
	#step 0: normalize
	#Find y_mean and y_std, use them to standardize data, and keep
	print("getting data...",flush=True) if verbose else 0
	p_mean = np.mean(qoi_samples)
	p_std = np.std(qoi_samples)
	yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_samples]
	d_mean = np.mean(y_samples, axis=0)
	d_std = np.std(y_samples, axis=0)
	yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_samples)]
	#print(p_mean)
	#print(p_std)
	#print(d_mean)
	#print(d_std)
	y_mean = [p_mean]+list(d_mean)
	y_std = [p_std]+list(d_std) #TODO put these into all of the GaussianMixtures

	p_dimension = 1 #len(yp_sample[0]) #1
	d_dimension = len(yd_sample[0]) #3

	data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_samples)])

	#step 1: train 
	#Get a a Gaussian mixture model from the push-forward of the prior through yd and yp

	#need to use a "Bayesian information criterion approach [8]"
	#or, a "number of components to balance the maximization of likelihood of data and the Bayesian information criterion [4]"
	#which one??? the latter seems fine
	#ok, generally to do this, we'll search bottom-up, increasing n_comp until we find a minimal BIC
	gmm = None
	if ncomp == 0:
		print("calculating ideal ncomp...",flush=True) if verbose else 0
		curr_gmm = GaussianMixtureNormalized(n_components=ncomp_start, max_iter=1000, n_init=5, reg_covar=1e-5, standardized_mean=y_mean, standardized_std=y_std) if careful else GaussianMixtureNormalized(n_components=ncomp_start, standardized_mean=y_mean, standardized_std=y_std)
		curr_gmm.fit(data)
		curr_bic = curr_gmm.bic(data)
		print(curr_bic, flush=True) if verbose==2 else 0
		next_gmm = GaussianMixtureNormalized(n_components=ncomp_start+1, max_iter=1000, n_init=5, reg_covar=1e-5, standardized_mean=y_mean, standardized_std=y_std) if careful else GaussianMixtureNormalized(n_components=ncomp_start+1, standardized_mean=y_mean, standardized_std=y_std)
		next_gmm.fit(data)
		next_bic = next_gmm.bic(data)
		print(next_bic, flush=True) if verbose==2 else 0
		ncomp = 1
		while next_bic < curr_bic:
			curr_bic = next_bic
			curr_gmm = next_gmm
			ncomp += 1
			if careful:
				try:
					next_gmm = GaussianMixtureNormalized(n_components=ncomp+1, max_iter=1000, n_init=5, reg_covar=1e-5, standardized_mean=y_mean, standardized_std=y_std).fit(data)
				except:
					#Exceptions are likely to happen if the number of components is larger than the sample size can support without introducing singularities in the covariance
					next_bic = curr_bic * 2 #flunk out to end the loop if it fails to converge
					print("(GMM with",ncomp+1,"components failed to converge)") if verbose else 0
					continue
			else:
				next_gmm = GaussianMixtureNormalized(n_components=ncomp+1, standardized_mean=y_mean, standardized_std=y_std).fit(data)
			next_bic = next_gmm.bic(data)
			print(next_bic, flush=True) if verbose else 0

		print("Number of components that minimizes BIC is",ncomp, flush=True) if verbose else 0
		gmm = curr_gmm
		
	else:
		#Allow someone to override this process if they think they know how many components they want
		if careful:
			gmm = GaussianMixtureNormalized(n_components=ncomp, max_iter=1000, n_init=5, reg_covar=1e-5, standardized_mean=y_mean, standardized_std=y_std).fit(data)
		else:
			gmm = GaussianMixtureNormalized(n_components=ncomp, standardized_mean=y_mean, standardized_std=y_std).fit(data)
		print("Using provided number of components,",ncomp,flush=True) if verbose else 0
		
	return gmm


#as above, but take in a warm-started gmm and continue to train it
def gbi_train_model_warm(qoi_samples, y_samples, gmm, verbose=0, careful=False):
	ncomp = len(gmm.weights_)
	gmm.warm_start = True
	
	print("getting data...",flush=True) if verbose else 0
	p_mean = np.mean(qoi_samples)
	p_std = np.std(qoi_samples)
	yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_samples]
	d_mean = np.mean(y_samples, axis=0)
	d_std = np.std(y_samples, axis=0)
	yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_samples)]
	#print(p_mean)
	#print(p_std)
	#print(d_mean)
	#print(d_std)
	y_mean = [p_mean]+list(d_mean)
	y_std = [p_std]+list(d_std) #TODO put these into all of the GaussianMixtures

	p_dimension = 1 #len(yp_sample[0]) #1
	d_dimension = len(yd_sample[0]) #3

	data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_samples)])

	#step 1: train 
	#Get a a Gaussian mixture model from the push-forward of the prior through yd and yp
	if careful:
		GaussianMixtureNormalized(n_components=ncomp+1, max_iter=1000, n_init=5, reg_covar=1e-5, standardized_mean=y_mean, standardized_std=y_std).fit(data)
	else:
		gmm = GaussianMixtureNormalized(n_components=ncomp, standardized_mean=y_mean, standardized_std=y_std).fit(data)
	print("Using provided number of components,",ncomp,flush=True) if verbose else 0
	
	return gmm


def gbi_precalc_Sigdd(gmm, p_dim=1):
	ncomp = len(gmm.weights_)
	#just invert the Sig_dd matrices once, saves a lot of time!
	Sig = [np.array(cov_k) for cov_k in gmm.covariances_]
	Sig_dd = [Sig_k[p_dim:, p_dim:] for Sig_k in Sig]
	inv_Sig_dd = [np.linalg.inv(Sig_dd_k) for Sig_dd_k in Sig_dd]
	
	#because the det is so small, I'm calling for the log of the det directly. slogdet requires careful unpacking
	signs_logdets = [np.linalg.slogdet(Sig_dd_k) for Sig_dd_k in Sig_dd] #returns (sign, logdet)
	logdet_Sig_dd = [None]*ncomp
	for k in range(ncomp):
		sign, logdet = signs_logdets[k]
		if sign == 1:
			logdet_Sig_dd[k] = logdet.item()
		elif sign < 0:
			print("Error: Sig_dd apparently wasn't positive semi-definite?")
			sys.exit()
		else:
			print("Error: Sig_dd has a 0 determinant, that's unexpected")
			sys.exit()
	
	return inv_Sig_dd, logdet_Sig_dd
	
#using decimal instead of mpmath
def gbi_condition_model(gmm, Yd_raw, inv_Sig_dd_precalc=None, logdet_Sig_dd_precalc=None, verbose=0):
	ncomp = gmm.n_components
	#step 2
	#Now we have our data, and we find the posterior predictive from that
	print("calculating posterior predictive...",flush=True) if verbose else 0
	
	##NOTE that we need to normalize this input data the same way that GaussianMixtureNormalized is:
	yd_means = gmm.standardized_mean[1:]
	yd_stds = gmm.standardized_std[1:]
	Yd_standard = [(y-yd_means[j])/yd_stds[j] for j,y in enumerate(Yd_raw)]
	Yd = np.array(Yd_standard)
	
	#get key parameters from the GMM
	mu = [np.array(mu_k) for mu_k in gmm.means_]
	Sig = [np.array(cov_k) for cov_k in gmm.covariances_]
	alpha = np.array(gmm.weights_)
	p_dimension = len(mu[0]) - len(Yd) #scalar, this is usually 1
	ymean_p = [mu_k[:p_dimension] for mu_k in mu] #column
	ymean_d = [mu_k[p_dimension:] for mu_k in mu] #column (kinda)
	Sig_pp = [Sig_k[:p_dimension, :p_dimension] for Sig_k in Sig]
	Sig_pd = [Sig_k[:p_dimension, p_dimension:] for Sig_k in Sig]
	Sig_dp = [Sig_k[p_dimension:, :p_dimension] for Sig_k in Sig]
	Sig_dd = [Sig_k[p_dimension:, p_dimension:] for Sig_k in Sig]
	
	###These things can be precalculated, to save time. check if they have been
	if inv_Sig_dd_precalc and logdet_Sig_dd_precalc:
		inv_Sig_dd = inv_Sig_dd_precalc
		logdet_Sig_dd = logdet_Sig_dd_precalc
	else: #no half-measures!
		inv_Sig_dd, logdet_Sig_dd = gbi_precalc_Sigdd(gmm, p_dim=1)

	#parameters for the new GMM:
	column_diff = [np.array(Yd - ymean_d[k]).reshape(len(Yd), 1) for k in range(ncomp)]
	exponent = [-0.5 * (column_diff[k].T @ inv_Sig_dd[k] @ column_diff[k]).item() for k in range(ncomp)] #TODO matrix mult isnt workin here
	logB1 = [(-ncomp/2.0)*np.log(alpha[k]*2*math.pi) 
			- 0.5 * logdet_Sig_dd[k]
			+ exponent[k]
			for k in range(ncomp)]
	#B1 = [np.exp(logbk) for logbk in logB1] #this is broken; would need to use Decimal to do this without loss of precision
	#B0 = math.fsum(B1) #this lets us sum floating-point numbers without loss of precision
	
	#scipy.special.logsumexp does log( SUMi exp(ai) )
	#-- similar to np.logaddexp, which does log(exp(a) + exp(b) safely
	logB0 = logsumexp(logB1)
	#lnB0 = ln [SUMk exp(lnBk)], equivalent to B0 = sum(Bk)
	
	#because we're dividing really small numbers, sometimes we get 1e-10 terms (or 1e-100 !), but those aren't harmful so leave em in
	#beta = [B1[k] / B0 for k in range(ncomp)] 
	beta = [np.exp(logB1[k] - logB0) for k in range(ncomp)]
	mu_Yd = [ymean_p[k] + Sig_pd[k] @ inv_Sig_dd[k] @ column_diff[k] for k in range(ncomp)]
	Sig_Yd = [Sig_pp[k] - Sig_pd[k] @ inv_Sig_dd[k] @ Sig_dp[k] for k in range(ncomp)]
	
	#convert back to np arrays, annoyingly: #TODO is this necessary?
	beta = [np.array(x) for x in beta]
	mu_Yd = [np.array(x) for x in mu_Yd]
	Sig_Yd = [np.array(x) for x in Sig_Yd]
	
	##LASTLY its very important that we de-normalize this input data according to the normalization of Q in the GaussianMixtureNormalized
	#should just involve multiplying the mean and std
	yp_mean = gmm.standardized_mean[0]
	yp_std = gmm.standardized_std[0]
	mu_Yd = [yp_std*mu + yp_mean for mu in mu_Yd] #xnorm = (x-xmean)/xstd
	Sig_Yd = [(yp_std**2)*var for var in Sig_Yd] #rescaling stddev=1 back to yp_std
	
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
	
def gbi_var_of_conditional_pp(gmm, Yd, verbose=0, inv_Sig_dd_precalc=None, logdet_Sig_dd_precalc=None):
	beta, mu_Yd, Sig_Yd = gbi_condition_model(gmm, Yd, inv_Sig_dd_precalc=inv_Sig_dd_precalc, logdet_Sig_dd_precalc=logdet_Sig_dd_precalc, verbose=verbose)
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

"""
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
"""