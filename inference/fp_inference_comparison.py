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
from problems.fp_verification.fp_vv_opt import *
from inference.lieberman_willcox_2014_sketch import *
from obed.mcmc import *
from uq.uncertainty_propagation import *

###Now try the focal plane case:
d_historical = [
				20,   #t_gain
				30,   #I_gain
				1,	#n_meas_rn
				8,	#d_num
				9600, #d_max
				2	 #d_pow   #approx
			   ]


#####Do the GMM posterior inference
def expt_fp(theta):
	yd = fp.eta(theta, d_historical)
	return yd

def pred_fp(theta):
	yp = fp.H(theta)
	return yp
	
def theta_prior_fp(num):
	return fp.prior_rvs(num)
	
theta_nominal = [1.1, 2.5, .001]
y_nominal = fp_likelihood_fn(dict(zip(fp.theta_names, theta_nominal)), dict(zip(fp.d_names, d_historical)), dict(zip(fp.x_names, fp.x_default)), err=False)
Yd = y_nominal
print("Simulated data:", Yd, flush=True)
print("QoI prediction given that data assuming no error:", pred_fp(theta_nominal), flush=True)

a,b,c = inference_predictive_posterior(expt_fp, pred_fp, theta_prior_fp, Yd, nsamples=10000, verbose=2)

#####Do the MCMC posterior, then the QoI push-through
prop_width = [0.007931095589546992, 0.018919515987306634, 0.00017949891623054683]
mcmc_trace,_,_ = mcmc_multigauss_likelihood(y_nominal, d_historical, proposal_fn_norm, prop_width, fp.eta, fp.prior_rvs, fp.prior_pdf_unnorm, n_mcmc=6000, n_pde=1000, burnin=300, lag=1, doPlot=True, legend=fp.theta_names, doPrint=True)
H_posterior = [fp.H(tt) for tt in mcmc_trace]


#plot both
plot_predictive_posterior(a, b, c, 0, 7, drawplot=False)
uncertainty_prop_plot(H_posterior, c='orangered', xlab="Posterior QoI: Avg. Noise [e-]", rescaled=True)

