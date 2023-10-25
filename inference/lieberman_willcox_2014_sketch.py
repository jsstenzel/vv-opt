#trying to implement the ideas in Lieberman & Willcos 2014

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

sys.path.append('..')
#focal plane
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

samples = theta_prior_sample(10000)
yd_sample = [expt_process(theta) for theta in samples]
yp_sample = [pred_process(theta) for theta in samples]

#actually i think ill use sklearn for this
data = np.array([np.hstack([yd_sample[i],yp_sample[i]]) for i,_ in enumerate(samples)])

#need to use a "Bayesian information criterion approach [8]"
#or, a "number of components to balance the maximization of likelihood of data and the Bayesian information criterion [4]"
#which one??? the latter seems fine
#ok, generally to do this, we'll search bottom-up, increasing n_comp until we find a minimal BIC
for i in [1,2,3,4,5]:
	ncomp = i
	gmm = GaussianMixture(n_components=ncomp).fit(data)
	bic = gmm.bic(data)
	print(bic)


