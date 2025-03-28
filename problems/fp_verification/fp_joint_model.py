import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import csv
import joblib

sys.path.append('../..')
#focal plane
from problems.fp_verification.fp_problem import *
#analysis
from obed.obed_gbi import *
from uq.uncertainty_propagation import *
from uq.plotmatrix import *
#from uq.sensitivity_analysis import *
from opt.ngsa import *

################################
#Useful definitions
################################

req = 4.38 #max noise

theta_nominal = [1.1, 2.5, .001]
QoI_nominal = fp.H(theta_nominal)

d_historical = [
				20,   #t_gain
				30,   #I_gain
				1,	#n_meas_rn
				8,	#d_num
				9600, #d_max
				2	 #d_pow   #approx
			   ]
			   
d_best = [
				600,   #t_gain
				100,   #I_gain
				50,	#n_meas_rn
				100,	#d_num
				12000, #d_max
				3	 #d_pow   #approx
		]
			   
d_worst = [
				1,   #t_gain
				1,   #I_gain
				1,	#n_meas_rn
				2,	#d_num
				1, #d_max
				0.1	 #d_pow   #approx
		]
		
y_nominal = fp_likelihood_fn(dict(zip(fp.theta_names, theta_nominal)), dict(zip(fp.d_names, d_historical)), dict(zip(fp.x_names, fp.x_default)), err=False)
#print(y_nominal)

################################
#Analysis functions
################################

def fp_vv_nominal():
	print("QoI requirement:", req)
	print("Nominal QoI:", QoI_nominal)


###uncertainty analysis
def fp_vv_UP_QoI(req):
	#uncertainty propagation of HLVA
	uq_thetas = fp.prior_rvs(10000)
	Qs = [fp.H(theta) for theta in uq_thetas]
	uncertainty_prop_plot([theta[0] for theta in uq_thetas], xlab="Gain [ADU/e-]")
	uncertainty_prop_plot([theta[1] for theta in uq_thetas], xlab="Read noise [e-]")
	uncertainty_prop_plot([theta[2] for theta in uq_thetas], xlab="Dark current [e-/s]")
	uncertainty_prop_plot(Qs, xlab="QoI: Avg. Noise [e-]", vline=[req])

	#prob of meeting req along priors:
	count_meetreq = 0
	for Q in Qs:
		if Q <= req:
			count_meetreq += 1
	prob_meetreq = count_meetreq / len(Qs)
	print("Probability of meeting requirement given priors:", prob_meetreq)

"""
#sensitivity analysis of HLVA
def fp_vv_SA_QoI():
	#it'll be straightforward to see the dependence of QoI on theta
	Si = sobol_saltelli(fp.H, 
						2**5, #SALib wants powers of 2 for convergence
						var_names=fp.theta_names, 
						var_dists=[prior[0] for prior in fp.priors], 
						var_bounds=[prior[1] for prior in fp.priors], 
						conf = 0.95, doSijCalc=False, doPlot=True, doPrint=True)
"""

#Uncertainty analysis of the experiment models
def fp_vv_UP_exp(dd, savefig=False):
	print("Likelihood distribution for nominal historical case:",flush=True)
	tt = fp.prior_rvs(1); print(tt)
	ysample_nominal = [fp.eta(tt, dd) for _ in range(10000)]
	uncertainty_prop_plots(ysample_nominal, xlabs=["Y0","Y1","Y2"], saveFig='UP_exp_nominal' if savefig else '')
	likelihood,_ = general_likelihood_kernel(ysample_nominal)
	#kde_plot(likelihood, ysample_nominal, plotStyle='together') #needs fixing?
	
	#Also plot the y's generated by the joint distribution p(y|theta,d)p(theta)
	print("Experiments simulated from the joint distribution:",flush=True)
	uq_thetas = fp.prior_rvs(10000)
	uq_ys = [fp.eta(theta, dd) for theta in uq_thetas]
	uncertainty_prop_plot([y[0] for y in uq_ys], xlab="Y0 (joint distribution)", c='orchid', saveFig='UP_joint_y0' if savefig else '')
	uncertainty_prop_plot([y[1] for y in uq_ys], xlab="Y1 (joint distribution)", c='orchid', saveFig='UP_joint_y1.png' if savefig else '')
	uncertainty_prop_plot([y[2] for y in uq_ys], xlab="Y2 (joint distribution)", c='orchid', saveFig='UP_joint_y2' if savefig else '')


###obed analysis -- naive multiprocessing approach
def fp_vv_obed_1machine():
	print("Generating kernel",flush=True)
	n_kde_axis = 47 #47^3 equals about 10^5
	kde_gains = np.linspace(0,3,n_kde_axis)
	kde_rn = np.linspace(1,4,n_kde_axis)
	kde_dc = np.linspace(0,.01,n_kde_axis)
	kde_thetas = np.vstack((np.meshgrid(kde_gains, kde_rn, kde_dc))).reshape(3,-1).T #thanks https://stackoverflow.com/questions/18253210/creating-a-numpy-array-of-3d-coordinates-from-three-1d-arrays
	kde_ys = [fp.eta(theta, d_historical) for theta in kde_thetas]
	likelihood_kernel, kde_ythetas = general_likelihood_kernel(kde_thetas, kde_ys)

	print("starting obed",flush=True)
	prop_width = [0.007167594573520732, 0.17849464019335232, 0.0006344271319903282] #stddev from last study

	U, U_list = U_probreq_multi(d_historical, fp, proposal_fn_norm, prop_width, likelihood_kernel, req=3.0, n_mc=1000, n_mcmc=2000, burnin=300, lag=1, doPrint=True)
	print(U)
	print(U_list)
	uncertainty_prop_plot(U_list, c='royalblue', xlab="specific U", saveFig='OBEDresult')
	
###obed analysis -- robust for sbatch, one step of the monte carlo
def fp_vv_obed_cluster(d, kernel_pkl):
	print("starting obed",flush=True)
	prop_width = [0.007167594573520732, 0.17849464019335232, 0.0006344271319903282] #stddev from last study

	U = U_probreq_1step(d_historical, fp, proposal_fn_norm, prop_width, kernel_pkl, req=3.0, n_mcmc=2000, burnin=300, lag=1, doPrint=True)

def fp_vv_plot_obed_results(file):
	u = []
	H_mean = []
	H_stddev = []
	with open(file) as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			u.append(float(row[0]))
			H_mean.append(float(row[4]))
			H_stddev.append(float(row[5]))
	
	uncertainty_prop_plot(u, c='royalblue', xlab="specific U from MC")
	uncertainty_prop_plot(H_mean, c='royalblue', xlab="H means from MC")
	uncertainty_prop_plot(H_stddev, c='royalblue', xlab="H stddevs from MC")
	print("MC-mean probability of meeting requirement:",np.mean(u))
	print("MC-mean of H_posterior mean:",np.mean(H_mean))
	print("   Probability of MC-mean meeting requirement:",len([1 for hmean in H_mean if hmean<3.0])/len(H_mean))
	print("MC-mean of H_posterior stddev:",np.mean(H_stddev))
	
	plt.scatter(H_mean, H_stddev)
	plt.show()
	
def fp_vv_gbi_test(d, Yd, N):		
	theta_train = fp.prior_rvs(N)
	qoi_train = [fp.H(theta) for theta in theta_train]
	y_train = [fp.eta(theta, d) for theta in theta_train]
	
	gmm = gbi_train_model(theta_train, qoi_train, y_train, verbose=2, ncomp=8)
	
	a,b,c = gbi_condition_model(gmm, Yd, verbose=2)
	
	plot_predictive_posterior(a, b, c, 0, 7, drawplot=True)
	
def fp_vv_gbi_rand_test(d, N):		
	theta_train = fp.prior_rvs(N)
	qoi_train = [fp.H(theta) for theta in theta_train]
	y_train = [fp.eta(theta, d) for theta in theta_train]
	
	gmm = gbi_train_model(theta_train, qoi_train, y_train, verbose=0, ncomp=8)
	
	truth_theta = fp.prior_rvs(1)
	measured_y = fp.eta(truth_theta, d)
	a,b,c = gbi_condition_model(gmm, measured_y, verbose=0)
	
	plot_predictive_posterior(a, b, c, 0, 7, drawplot=False, plotmean=True)
	plt.axvline(fp.H(truth_theta), c='blue')
	plt.show()

def fp_vv_obed_gbi(d):
	U, U_list = U_varH_gbi(d, fp, n_mc=10**5, n_gmm=10**5, doPrint=True)
	print(U)
	#print(U_list)
	uncertainty_prop_plot(U_list, c='royalblue', xlab="specific U")#, saveFig='OBEDresult')
	return U
	
def gbi_test(d):
	###get a validation set
	theta_val = problem.prior_rvs(N)
	d_val = [d for d in problem.sample_d(N)]
	qoi_val = [problem.H(theta) for theta in theta_val]
	y_val = [problem.eta(theta, d) for theta,d in zip(theta_val,d_val)]
	val_input = [yi+di for yi,di in zip(y_val, d_val)]
	val_output = np.array(qoi_val)

def fp_jm_test_plot(problem, N):
	###get a validation set
	theta_val = problem.prior_rvs(N)
	d_val = [d for d in problem.sample_d(N)]
	qoi_val = [problem.H(theta) for theta in theta_val]
	y_val = [problem.eta(theta, d) for theta,d in zip(theta_val,d_val)]
	val_input = [yi+di for yi,di in zip(y_val, d_val)]
	val_output = np.array(qoi_val)
	
	###plot the covariance?
	print(np.shape(val_input))
	print(np.shape(val_output))
	
	#plot just the y,d data
	#names = problem.y_names + problem.d_names
	#plotmatrix(np.array(val_input))
	
	#plot it alongside q
	names = problem.y_names + problem.d_names + ['q']
	plotmatrix(np.append(np.array(val_input),val_output[:,None], axis=1), names=names, rescale=False)
	
#i dont think my little test here is actually well-formed
def fp_jm_test_gmm(problem, N, N_val):
	###get a validation set
	theta_val = problem.prior_rvs(N_val)
	d_val = [d for d in problem.sample_d(N_val)]
	qoi_val = [problem.H(theta) for theta in theta_val]
	y_val = [problem.eta(theta, d) for theta,d in zip(theta_val,d_val)]
	val_input = [yi+di for yi,di in zip(y_val, d_val)]
	val_output = np.array(qoi_val)
	
	###get the model
	model = joint_model_gmm(problem, N, doPrint=True)
	
	###calculate what the model calculates for the validation set
	model_output = [gbi_sample_of_conditional_pp(model, vi) for vi in val_input]
	
	###compare
	validation_error = abs(model_output - val_output)
	print(np.mean(validation_error))
	uncertainty_prop_plot(validation_error, xlab="model prediction error", c='r')
	
#here i want to specifically compare how the new yd-q joint model performs compared to the y-q joint model from obed_gbi
#to do that at one data point, I need to pick a y and d (implicitly picking a theta)
#then, I need to train the ydq joint model with gbi_train_model(joint_output, joint_output, joint_input, ncomp=0, verbose=doPrint)
#and I need to train the yq joint model with gbi_train_model(theta_train, qoi_train, y_train, ncomp=ncomp)
def fp_jm_test_gmm_var(problem, N_yq, N_ydq, N_val, valplots=True, d_fixed=[]):
	###get a training set
	theta_train = problem.prior_rvs(N_ydq)
	d_train = [d for d in problem.sample_d(N_ydq)]
	qoi_train = [problem.H(theta) for theta in theta_train]
	y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
	train_input = [yi+di for yi,di in zip(y_train, d_train)]
	train_output = np.array(qoi_train)
	
	ydq_joint_gmm = gbi_train_model(train_output, train_output, train_input, ncomp=0, verbose=True)
	
	###get a different training set
	#theta_train = problem.prior_rvs(N_yq)
	#d_train = [d for d in problem.sample_d(N_yq)]
	#qoi_train = [problem.H(theta) for theta in theta_train]
	#y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
	#train_input = [yi+di for yi,di in zip(y_train, d_train)]
	#train_output = np.array(qoi_train)
	
	###get the yq joint at d
	if len(d_fixed)==0:
		d_fixed = problem.sample_d(1)
		print(d_fixed)
	
	y_train_d_fixed = [problem.eta(theta, d_fixed) for theta in theta_train]
	yq_joint_gmm = gbi_train_model(theta_train, qoi_train, y_train_d_fixed, ncomp=0, verbose=True)
	
	###then compare posterior predictive variance of both models at d
	theta_val = problem.prior_rvs(N_val)
	qoi_val = [problem.H(theta) for theta in theta_val]
	y_val = [problem.eta(theta, d_fixed) for theta in theta_val]
	val_input = [yi+d_fixed for yi in y_val]
	
	yq_model_eval = []
	ydq_model_eval = []
	for yi,vi in zip(y_val, val_input):
		#Here, I'm essentially using two different models to calculate the mean varH over theta and y(theta,d_fixed)
		#or in other words, I'm using two different models to calculate U(d) where n=N_val
		beta_y, mu_y, Sig_y = gbi_condition_model(yq_joint_gmm, yi)
		beta_v, mu_v, Sig_v = gbi_condition_model(ydq_joint_gmm, vi)
		
		if valplots:
			plot_predictive_posterior(beta_y, mu_y, Sig_y, 0, 8, compplot=False, drawplot=False, plotmean=True)
			plot_predictive_posterior(beta_v, mu_v, Sig_v, 0, 8, compplot=False, drawplot=True, plotmean=True, maincolor='darkgreen')

		yq_model_eval.append(gbi_gmm_variance(beta_y, mu_y, Sig_y))
		ydq_model_eval.append(gbi_gmm_variance(beta_v, mu_v, Sig_v))
	
	###compare
	if valplots:
		model_mean_diff = [abs(u1 - u2) for u1,u2 in zip(yq_model_eval, ydq_model_eval)]
		uncertainty_prop_plot(model_mean_diff, xlab="model prediction error", c='r')
		print("U diff:",np.mean(model_mean_diff))
	return np.mean(yq_model_eval), np.mean(ydq_model_eval)
	
#make a new version of the above that calculates this at many d
def fp_jm_test_gmm_U(problem, N_train, N_mc, N_opt):
	###get a training set
	theta_train = problem.prior_rvs(N_train)
	d_train = [d for d in problem.sample_d(N_train)]
	qoi_train = [problem.H(theta) for theta in theta_train]
	y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
	train_input = [yi+di for yi,di in zip(y_train, d_train)]
	train_output = np.array(qoi_train)
	
	ydq_joint_gmm = gbi_train_model(train_output, train_output, train_input, ncomp=0, verbose=False)
	
	###get the yq joint at d
	d_val = [d for d in problem.sample_d(N_opt)]
	val_err = []
	for d_fixed in d_val:
		print(d_fixed)
		y_train_d_fixed = [problem.eta(theta, d_fixed) for theta in theta_train]
		yq_joint_gmm = gbi_train_model(theta_train, qoi_train, y_train_d_fixed, ncomp=0, verbose=False)
		
		###then evaluate both models at d
		theta_val = problem.prior_rvs(N_mc)
		qoi_val = [problem.H(theta) for theta in theta_val]
		y_val = [problem.eta(theta, d_fixed) for theta in theta_val]
		val_input = [yi+d_fixed for yi in y_val]
		
		#Here, I'm essentially using two different models to calculate the mean varH over theta and y(theta,d_fixed)
		#or in other words, I'm using two different models to calculate U(d) where n=N_mc
		U_yq = np.mean(np.array([gbi_var_of_conditional_pp(yq_joint_gmm, yi) for yi in y_val]))
		U_ydq = np.mean(np.array([gbi_var_of_conditional_pp(ydq_joint_gmm, vi) for vi in val_input]))
		
		###compare
		validation_error = abs(U_yq - U_ydq)/U_yq
		val_err.append(validation_error)
	
	print(val_err, flush=True)
	uncertainty_prop_plot(val_err, xlab="model prediction error", c='r')
	
def fp_vv_obed_gbi_joint(d, gmm, n):
	U, U_list = U_varH_gbi_joint(d, fp, gmm, n_mc=n, doPrint=True)
	print(U)
	print(U_list)
	uncertainty_prop_plot(U_list, c='royalblue', xlab="specific U")#, saveFig='OBEDresult')
	return U


if __name__ == '__main__':  
	#fp_vv_nominal()
	
	#fp_vv_UP_QoI(4.38)
	
	#fp_vv_SA_QoI()
	
	#fp_vv_UP_exp(d_historical)
	#fp_vv_UP_exp(d_best, False)
	#fp_vv_UP_exp(d_worst, False)
	
	#fp_vv_gbi_test(d_historical, y_nominal, 10**6)
	
	#fp_vv_gbi_rand_test(d_historical, 10**4)
	
	#U_hist = fp_vv_obed_gbi(d_historical)
		
	#fp_jm_test_plot(fp, 10**3)
	
	#fp_jm_test_gmm(fp, 50, 10**3)
	#fp_jm_test_gmm(fp, 10**2, 10**3)
	#fp_jm_test_gmm(fp, 10**3, 10**3)
	
	#fp_jm_test_gmm_var(fp, 10**3, 10**2)
	#fp_jm_test_gmm_var(fp, 10**4, 10**2)
	#fp_jm_test_gmm_var(fp, 10**5, 10**2)
	#this shows performance getting worse for higher N
	#I think thats maybe not surprising, both GMMs are getting better, i guess at different rates
	#I could freeze N for the yq one, but even then that may not be terribly informative, what if the ydq one exceeds it?
	#fp_jm_test_gmm_var(fp, N_yq=10**3, N_ydq=4*10**3, N_val=10)
	
	#convergence plot
	nlist = [10**3, 4*10**3, 10**4, 4*10**4, 10**5, 4*10**5, 10**6, 4*10**6]
	d_fixed = fp.sample_d(1)
	ulist_yq = []
	ulist_ydq = []
	for n in nlist:
		u_yq, u_ydq = fp_jm_test_gmm_var(fp, N_yq=n, N_ydq=n, N_val=10**3, d_fixed=d_fixed, valplots=False)
		ulist_yq.append(u_yq)
		ulist_ydq.append(u_ydq)
	plt.plot(nlist, ulist_yq, c='k')
	plt.plot(nlist, ulist_ydq, c='darkgreen')
	plt.show()
	
	#What is actually a good test here?
	#Hard to test against "truth" because that involves MCMC, which is not my friend anymore
	#Nannapaneni 2020 suggests comparing pareto fronts
	#ah fuck worry about it later then!
	
	#fp_jm_test_gmm_U(fp, N_train=10**3, N_mc=10**5, N_opt=25)
	
	#train the GMM of our BN
	#ydq_joint_gmm = joint_model_gmm(fp, 10**6, doPrint=True, doPlot=False)
	#actually, i should do this totally separately, so i can manually record the properties of the gmm for posterity & repetition
	# save
	#joblib.dump(ydq_joint_gmm, "gmm.pkl") 

	"""
	# load
	ydq_joint_gmm = joblib.load("gmm_10^6.pkl")
	
	#now use it
	fp_vv_obed_gbi_joint(d_historical, ydq_joint_gmm, 10**2)
	fp_vv_obed_gbi_joint(d_best, ydq_joint_gmm, 10**2)	
	fp_vv_obed_gbi_joint(d_worst, ydq_joint_gmm, 10**2)	
	"""
	
	"""
	#U_hist = fp_vv_obed_gbi(d_historical)
	U_hist = 0.4981609760099987
	C_hist = np.log(fp.G(d_historical))
	
	#U_best = fp_vv_obed_gbi(d_best)
	U_best = 0.48310825880372715
	C_best = np.log(fp.G(d_best))
	
	#U_worst = fp_vv_obed_gbi(d_worst)
	U_worst = 1.2584888138976271
	C_worst = np.log(fp.G(d_worst))

	pts = [[C_hist, U_hist, "d_hist"], [C_best, U_best, "d_max"], [C_worst, U_worst, "d_min"]]
	plot_ngsa2(costs, utilities, design_pts=pts, showPlot=True, savePlot=False, logPlotXY=[False,False])
	"""
