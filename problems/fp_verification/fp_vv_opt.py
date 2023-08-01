import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv

sys.path.append('../..')
#focal plane
from problems.problem_definition import *
from problems.fp_verification.fp_statistics import *
from problems.fp_verification.fp_experiment_models import *
#snslysis
from obed.obed_multivar import *
from uq.uncertainty_propagation import *
from uq.sensitivity_analysis import *

useQE = False

eta = fp_likelihood_fn
H = fp_hlva
Gamma = fp_cost_simple			   

#these priors are based on requirements that were met, see Camera Qual Report
theta_req_defs = [ 
					("gain", ["gamma_mv", [1.1,0.2**2]]), #mean, variance
					("rn",   ["gamma_mv", [2.5,0.25**2]]), 
					("dc",   ["gamma_mv", [0.001,.001**2]])
				]
#need to update with range somehow? These can't be negative

#for post-fabrication pre-testing priors, I should update the above priors
#with the Qual data from SN20001 and SN20003
#theta_qual_defs = [ 
#					("gain", ["gaussian", [1.0,0.2]]), #mean, stddev
#					("rn",   ["gaussian", [2.5,0.25]]),
#					("dc",   ["gaussian", [0.001,.002]]),
#				]

fp_y_defs = ["y_gain", "y_rn", "y_dc"]

fp_d_defs = [
				"t_gain", "I_gain",		   #gain
				"n_meas_rn",			   #rn
				"d_num", "d_max", "d_pow"  #dc
			 ]
	
_temp= -90+273.15 #K
_k = 1.380649e-23 #J / K
_c = 299792458 #m / s
_e0 = 22100 #eV
_m0 = 108.9049856 * 1.66054e-27 #cd109 atomic mass in kg
fp_x_defs = [
				#general
				("nx", 2048),
				("ny", 2048),
				("sigma_dc", .5), #e-/s #Estimate based on a consistency test performed on SN20006
				("mu_stray", 0), #e-/s #WAG for now
				("sigma_stray", .005), #WAG for now
				#gain
				("P_signal", 0.90), #Prob. of correctly identifying signal as event #WAG for now
				("P_noise", 0.01), #Prob. of incorrectly identifying noise/interference as event #WAG for now
				("T_ccd", _temp), #K
				("E0", _e0), #22.1 keV Cd-109 emission line
				("sigma_E", math.sqrt((_m0 * _c**2) / (_k*_temp*_e0**2))), #1/eV^2
				("w", 3.66 + 0.000615*(300-_temp)), #eV/e- #this has uncertainty. nuisance parameter?
				("activity_cd109", 5e-6), #Ci #radioactivity of sample
				("grade_size", 3), #3x3 event grade sizes
				("t_gain_setup", 1200), #WAG
				("t_gain_buffer", 5), #WAG
				#rn
				("t_rn", .1), #100 ms exposure
				("t_rn_buffer", 5), #WAG
				#dc
				("t_0", 0.1), #100ms baseline exposure assumed
				("t_dc_buffer", 5), #WAG
				#qoi
				("tau", 1800),
				#cost
				("testbed_setup", 1800), #WAG
				("C_engineer", 0.00694444444) #WAG $/s, from $25/hr
			]
			 
if useQE == True:
	eta = fp_qe_likelihood_fn
	H = fp_qe_hlva
	Gamma = fp_qe_cost			   

	#this prior is based on measurements for red CCD42-40
	theta_req_defs.append(
		("qe", ["funct_splines", [
					[(400,.25),(500,.45),(650,.75),(900,.45),(975,.05)], #measured data
					3, #order of interpolation
					.05, #y error on all points
					[350,975], #LLAMAS spectral range
					[0.0,1.0], #range of QE
			 ]]
		))

	fp_y_defs.append("y_qe")

	fp_d_defs.extend(["n_qe", "t_qe", "I_qe", "S_err"]) #qe design variables
		
	S_pd = Functional([(200,.12),(280,.1),(300,.125),(400,.185),(633,.33),(930,.5),(1000,.45),(1100,.15)]) #nm, A/W
	S_pd.set_xlim(350,975)
	S_pd.set_ylim(0,1)
	S_pd.spline_interp(3)
	fp_x_defs.extend([
					#general
					#gain
					#rn
					#dc
					#qe
					("S_pd", S_pd),   #representative photodiode, Hamamatsu S1337-1010BQ
					("S_pd_err", .01)  #mA/W
					#qoi
					#cost
					])


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
fp = ProblemDefinition(eta, H, Gamma, theta_req_defs, fp_y_defs, fp_d_defs, fp_x_defs)

if __name__ == '__main__':  
	###nominal case
	req = 3.0 #max noise
	print("QoI requirement:", req)

	theta_nominal = [1.1, 2.5, .001]
	QoI_nominal = fp.H(theta_nominal)
	print("Nominal QoI:", QoI_nominal)


	###uncertainty analysis
	if False:
		#uncertainty propagation
		uq_thetas = fp.prior_rvs(10000)
		Qs = [fp.H(theta) for theta in uq_thetas]
		uncertainty_prop_plot([theta[0] for theta in uq_thetas], xlab="Gain [ADU/e-]")
		uncertainty_prop_plot([theta[1] for theta in uq_thetas], xlab="Read noise [e-]")
		uncertainty_prop_plot([theta[2] for theta in uq_thetas], xlab="Dark current [e-/s]")
		uncertainty_prop_plot(Qs, xlab="QoI: Avg. Noise [e-]")

		#prob of meeting req along priors:
		count_meetreq = 0
		for Q in Qs:
			if Q <= req:
				count_meetreq += 1
		prob_meetreq = count_meetreq / len(Qs)
		print("Probability of meeting requirement given priors:", prob_meetreq)

	#sensitivity analysis
	if False:
		#it'll be straightforward to see the dependence of QoI on theta
		Si = sobol_saltelli(fp.H, 
							2**5, #SALib wants powers of 2 for convergence
							var_names=fp.theta_names, 
							var_dists=[prior[0] for prior in fp.priors], 
							var_bounds=[prior[1] for prior in fp.priors], 
							conf = 0.95, doSijCalc=False, doPlot=True, doPrint=True)


		#less strightforward to do it for the dependence of QoI on d
	
	d_historical = [
					20,   #t_gain
					30,   #I_gain
					1,	#n_meas_rn
					8,	#d_num
					9600, #d_max
					2	 #d_pow   #approx
				   ]
	
	#While we're here, let's look at model uncertainty of the experiment models
	if False:
		print("Likelihood distribution for nominal historical case:",flush=True)
		tt = fp.prior_rvs(1); print(tt)
		ysample_nominal = [fp.eta(tt, d_historical) for _ in range(10000)]
		uncertainty_prop_plot([y[0] for y in ysample_nominal], xlab="Y0")
		uncertainty_prop_plot([y[1] for y in ysample_nominal], xlab="Y1")
		uncertainty_prop_plot([y[2] for y in ysample_nominal], xlab="Y2")
		likelihood = general_likelihood_kernel(ysample_nominal)
		tt = fp.prior_rvs(1); print("likelihood", tt, likelihood(fp.eta(tt, d_historical)))
		tt = fp.prior_rvs(1); print("likelihood", tt, likelihood(fp.eta(tt, d_historical)))
		tt = fp.prior_rvs(1); print("likelihood", tt, likelihood(fp.eta(tt, d_historical)))
		tt = fp.prior_rvs(1); print("likelihood", tt, likelihood(fp.eta(tt, d_historical)))
		tt = fp.prior_rvs(1); print("likelihood", tt, likelihood(fp.eta(tt, d_historical)))
		likelihood_plot(likelihood, ysample_nominal)
		
		#Also plot the y's generated by the joint distribution p(y|theta,d)p(theta)
		print("Experiments simulated from the joint distribution:",flush=True)
		uq_thetas = fp.prior_rvs(10000)
		uq_ys = [fp.eta(theta, d_historical) for theta in uq_thetas]
		uncertainty_prop_plot([y[0] for y in uq_ys], xlab="Y0 (joint distribution)", c='orchid')
		uncertainty_prop_plot([y[1] for y in uq_ys], xlab="Y1 (joint distribution)", c='orchid')
		uncertainty_prop_plot([y[2] for y in uq_ys], xlab="Y2 (joint distribution)", c='orchid')

	###mcmc analysis
	y_nominal = fp_likelihood_fn(dict(zip(fp.theta_names, theta_nominal)), dict(zip(fp.d_names, d_historical)), dict(zip(fp.x_names, fp.x_default)), err=False)
	print(y_nominal)
	
	def proposal_fn_gamma(theta_curr):
		theta_prop = [0] * len(theta_curr)
		proposal_width = [.3**2,.3**2,.002**2]
		for i,_ in enumerate(theta_prop):
			#proposal dists are gammas, to match
			mean = abs(theta_curr[i])
			variance = proposal_width[i]
			alpha = mean**2 / variance
			beta = mean / variance
			theta_prop[i] = scipy.stats.gamma.rvs(size=1, a=alpha, scale=1.0/beta)[0]
		return theta_prop
		
	def proposal_fn_norm(theta_curr):
		theta_prop = [0] * len(theta_curr)
		proposal_width = [2,2,.01]
		for i,_ in enumerate(theta_prop):
			#proposal dists are gammas, to match
			mean = abs(theta_curr[i])
			stddev = proposal_width[i]
			theta_prop[i] = scipy.stats.norm.rvs(size=1, loc=mean, scale=stddev)[0]
		return theta_prop
		
	if True: #play with mcmc_multivar
		print("Generating kernel",flush=True)
		n_kde = 10**5
		theta_domains = [[0,3],[1,4],[0,.01]]
		kde_thetas = [[scipy.stats.uniform.rvs(size=1, loc=left, scale=right-left)[0] for left,right in theta_domains] for _ in range(n_kde)]
		kde_ys = [fp.eta(theta, d_historical) for theta in kde_thetas]
		likelihood_kernel = general_likelihood_kernel(kde_thetas, kde_ys)
		#likelihood_plot(likelihood_kernel, kde_ys) this guy is 6d, i need to generalize that algorithm first
		
		print("mcmc",flush=True)
		n_mcmc = 10**4
		mcmc_trace, arate, rrate = mcmc_kernel(y_nominal, likelihood_kernel, proposal_fn_norm, fp.prior_rvs, fp.prior_pdf_unnorm, n_mcmc, burnin=0, lag=1, doPlot=True, legend=fp.theta_names)
		print(arate, rrate)
	
	if False: #play with mcmc_nolikelihood
		print("mcmc",flush=True)
		n_kde = 10**4
		n_mcmc = 10**2
		mcmc_trace, arate, rrate = mcmc_nolikelihood(y_nominal, d_historical, fp.eta, proposal_fn_norm, fp.prior_rvs, fp.prior_pdf_unnorm, n_mcmc, n_kde, burnin=0, lag=1, doPlot=True, legend=fp.theta_names)
		print(arate, rrate)
	
	#obed analysis
	#U, U_list = U_probreq(d_historical, fp, maxreq=3.0, n_mc=100, n_mcmc=1000, burnin=0, lag=1)
	#print(U)
	#print(U_list)