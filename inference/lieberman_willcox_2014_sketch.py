#trying to implement the ideas in Lieberman & Willcos 2014

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

sys.path.append('..')
#focal plane
print("importing fp...",flush=True)
from problems.fp_verification.fp_problem import *

d_historical = [
				20,   #t_gain
				30,   #I_gain
				1,	#n_meas_rn
				8,	#d_num
				9600, #d_max
				2	 #d_pow   #approx
			   ]

def expt_process(theta):
	yd = fp.eta(theta, d_historical)
	return yd

def pred_process(theta):
	yp = fp.H(theta)
	return yp
	
def theta_prior_sample(num):
	return fp.prior_rvs(num)
	
def theta_prior_pdf(theta):
	return fp.prior_pdf_unnorm(theta)
	
#step 1: train 
#Get a a Gaussian mixture model from the push-forward of the prior through yd and yp

samples = theta_prior_sample(100000)
yp_sample = [pred_process(theta) for theta in samples]
yd_sample = [expt_process(theta) for theta in samples]

p_dimension = 1 #len(yp_sample[0]) #1
d_dimension = len(yd_sample[0]) #3

#actually i think ill use sklearn for this
print("getting data...",flush=True)
data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(samples)])

#need to use a "Bayesian information criterion approach [8]"
#or, a "number of components to balance the maximization of likelihood of data and the Bayesian information criterion [4]"
#which one??? the latter seems fine
#ok, generally to do this, we'll search bottom-up, increasing n_comp until we find a minimal BIC
"""
print("calculating idean ncomp...",flush=True)
for i in range(1,15):
	ncomp = i
	gmm = GaussianMixture(n_components=ncomp).fit(data)
	bic = gmm.bic(data)
	print(bic)
"""

ncomp = 7
gmm = GaussianMixture(n_components=ncomp).fit(data)

#step 2
#Now we have our data, and we find the posterior predictive from that
theta_nominal = [1.1, 2.5, .001]
y_nominal = fp_likelihood_fn(dict(zip(fp.theta_names, theta_nominal)), dict(zip(fp.d_names, d_historical)), dict(zip(fp.x_names, fp.x_default)), err=False)
Yd = y_nominal

#get key parameters from the GMM
mu = gmm.means_
Sig = gmm.covariances_
alpha = gmm.weights_
ymean_p = np.array([mu_k[:p_dimension] for mu_k in mu])
ymean_d = np.array([mu_k[p_dimension:] for mu_k in mu])
Sig_pp = np.array([Sig_k[:p_dimension, :p_dimension] for Sig_k in Sig])
Sig_pd = np.array([Sig_k[p_dimension:, :p_dimension] for Sig_k in Sig])
Sig_dp = np.array([Sig_k[:p_dimension, p_dimension:] for Sig_k in Sig])
Sig_dd = np.array([Sig_k[p_dimension:, p_dimension:] for Sig_k in Sig])

#parameters for the new GMM:
B1 = [alpha[k] * () * abs * np.exp()] #finish
B0 #finish
Beta = [B1[k] / B0[k] for k in range(ncomp)]
mu_Yd = [ymean_p[k] + Sig_pd[k] @ np.inv(Sig_dd[k]) @ (Yd - ymean_d[k]) for k in range(ncomp)]
Sig_Yd = [Sig_pp[k] - Sig_pd[k] @ np.inv(Sig_dd[k]) @ Sig_dp[k] for k in range(ncomp)]


#i need to validate this algorithm - check [17], original thesis










