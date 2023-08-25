import sys
import os
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('..')
from problems.problem_definition import *
from obed.mcmc import *
from uq.uncertainty_propagation import *
from obed.pdf_estimation import *

#This sketch is intended to test that the MCMC algorithm
#converges on the correct posterior distribution

#just using gaussians for this for now

PRIOR_MEAN = 0.0
PRIOR_STDDEV = 1.0

ETA_MEAN = 0.0
ETA_STDDEV = 2.0

def proposal_fn_norm(theta_curr, proposal_width):
	theta_prop = [0] * len(theta_curr)
	for i,_ in enumerate(theta_prop):
		#proposal dists are gammas, to match
		#mean = abs(theta_curr[i])
		mean = abs(theta_curr[i])
		stddev = proposal_width[i]
		theta_prop[i] = scipy.stats.norm.rvs(size=1, loc=mean, scale=stddev)[0]
	return theta_prop

def toy_eta(theta, d, x, err=True):
	t1 = theta["theta1"]
	random = scipy.stats.norm.rvs(loc=ETA_MEAN, scale=ETA_STDDEV, size=1)[0]
	return [t1 + random]
	
def toy_H(theta, x):
	return 0
	
def toy_Gamma(d, x):
	return 0

toy_theta_defs = [ 
					["theta1", ["gaussian", [PRIOR_MEAN, PRIOR_STDDEV]], "continuous"]
				 ]
toy_y_defs = ["y1"]
toy_d_defs = [
				["d1", ['uniform', [0.0, 0.0]], "continuous"]
			]
toy_x_defs = []
test = ProblemDefinition(toy_eta, toy_H, toy_Gamma, toy_theta_defs, toy_y_defs, toy_d_defs, toy_x_defs)

dd = [0] #doesnt matter here
yy = [PRIOR_MEAN + ETA_MEAN]

#####Check that my prior pdf looks good:
prior_probs = []
yis = np.arange(-3,3,.01)
for yi in yis:
	prior_probs.append(test.prior_pdf_unnorm([yi]))
plt.plot(yis, prior_probs)
plt.title("prior pdf (unnormalized)")
plt.show()

#####Check that the pdf estimation is working as expected:
tt = [PRIOR_MEAN]
log_probs = []
yis = np.arange(-8,8,.01)
for yi in yis:
	log_probs.append(eta_multigaussian_logpdf([yi], tt, dd, test.eta, 3000))
plt.plot(yis, [np.exp(lp) for lp in log_probs])
plt.title("posterior pdf estimate")
plt.show()


#####Then, run the MCMC
prop_width = [2/np.sqrt(5)]
#mcmc_trace,acceptance_rate,_ = mcmc_multigauss_likelihood(yy, dd, proposal_fn_norm, prop_width, toy.eta, toy.prior_rvs, toy.prior_pdf_unnorm, n_mcmc=20, n_pde=1000, burnin=0, lag=1, doPlot=True, legend=toy.theta_names, doPrint=True)

mcmc_trace,acceptance_rate,_ = mcmc_multigauss_likelihood(yy, dd, proposal_fn_norm, prop_width, test.eta, test.prior_rvs, test.prior_pdf_unnorm, n_mcmc=10000, n_pde=1000, burnin=0, lag=1, doPlot=True, legend=test.theta_names, doPrint=True)
print("Acceptance rate:",acceptance_rate)

print("mean, stddev, covariance of posterior sample:")
means, stddevs, cov = mcmc_analyze(mcmc_trace,doPlot=True)
print(means)
print(stddevs)
print(cov)
uncertainty_prop_plots(mcmc_trace, c='limegreen', xlabs=["posterior"])