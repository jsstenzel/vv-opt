#Building on L&W 2014 and goal_based_inference.py
#Here, I am creating methods for easy and efficient offline training of a GMM(y, d, Q)
#And the methods to sample iteratively, save and load samples, and show convergence of training
#The idea is that, at the end, I can create a GMM definition that itself can be saved and loaded (a list of means, vars, and coefficients)
#Compare to the general approach in saltelli_gsa

import os
import sys
import math
import csv
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('..')
from inference.bn_modeling import *
#from inference.goal_based_inference import *
from uq.plotmatrix import *
from uq.uncertainty_propagation import *

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
###MODEL EXAMINATION

bn_measure_model_mse -- a flawed algorithm for calculating MSE for one GMM. Revisit when I have MCMC working for real validation

bn_measure_stability_convergence -- a flawed algorithm trying to look at convergence of Q-inference on a validation y set

bn_compare_model_covariance -- a nice visual inspection to see how similar data and model covariances are

bn_evaluate_model_likelihood -- find avg loglikelihood or score of a single GMM

###NAIVE CONVERGENCE TESTING

bn_measure_likelihood_convergence -- train/save/load a series of N-GMMs, and measure convergence on the training data

bn_measure_validation_convergence -- does the same, but on validation data

###N_COMP EVALUATION

bn_train_evaluate_ncomp -- train/save/load a series of ncomp-GMMs, print the BIC

bn_train_evaluate_ncomp_plot -- dirty function to manually plot an ncomp-series of BICs

bn_train_evaluate_ncomp_sanitycheck -- similar to bn_train_evaluate_ncomp, loads a series of GMMs, plots BIC and likelihood and penalty terms

###BOOTSTRAPPED CONVERGENCE TESTING

bn_train_convergence_confidence -- train and save a series of bootstrapped GMMs with ncomp and Ntrain

bn_evaluate_convergence_confidence -- load that series, evaluate confidence interval of the validation data score
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	
#This calculates MSE between the truth and model
#Note that of course, conditioning the GMM gets you the posterior predictive of Q
#So what I actually want to compare here is how well Var[GMM(y,d)] compares to Var[Q]
#NOTE that this doesn't compare Var[GMM(y,d)] to Var[p(Q|y)], to do that I need to do inference, such as with MCMC
#This is a little busted; we're comparing Q data prior to Q model posterior which is apples and pears
#TODO someday implement MCMC in here to compare data MCMC posterior to model GBI posterior
def bn_measure_model_mse(problem, gmm, N, doPrint=True):
	N_val = N if N>0 else 1000
	if doPrint:
		print("Drawing",N_val,"validation samples...",flush=True)
	
	###Draw N_val totally new validation samples
	theta_samples = problem.prior_rvs(N_val)
	d_samples = problem.sample_d(N_val)
	
	###################################
	if doPrint:
		print("Evaluating truth values...",flush=True)
	###Evaluate the validation set for "truth" case
	#for each theta,d get a bunch of subsamples of Q, and calculate the variance of that
	subsample_N = N_val #i suppose?
	Qvar_true = []
	for i in range(N_val):
		qoi_train = [problem.H(theta_samples[i]) for i in range(subsample_N)]
		#y_train = [problem.eta(theta_samples[i], d_samples[i]) for _ in range(subsample_N)]
		Qvar_true.append(np.var(qoi_train, ddof=0))  #not ddof=1, since I want to compare with the conditioned GMM, which uses the "actual" variance!
	
	###################################
	if doPrint:
		print("Evaluating model values...",flush=True)
	###Evaluate the validation set with the GMM
	y_samples = [problem.eta(theta_samples[j], d_samples[j]) for j in range(N_val)]
	#this introduces uncertainty that is not present in the truth case. To calculate real MSE, i would need to actually calculate posterior
		
	v = [y+d for y,d in zip(y_samples, d_samples)] #turn into array
	Qvar_gmm = [gbi_var_of_conditional_pp(gmm, vi, verbose=2) for vi in v] #turn into array
	
	SE = [(truth - model)**2 for truth,model in zip(Qvar_true, Qvar_gmm)]
	MSE = np.mean(SE)
	if doPrint:
		print("MSE for",N_val,"validation samples:",MSE, flush=True)
		
	return MSE


def bn_measure_stability_convergence(problem, big_savefile, N_val, doPrint=True):
	N_val = N_val if N_val>0 else 5		
	N_list = [1000,4000,10000,40000,100000,400000]
	
	###Draw N_val totally new validation samples
	#TODO i probably should not be reusing validation samples.
	if doPrint:
		print("Drawing",N_val,"validation samples...",flush=True)
	theta_samples = problem.prior_rvs(N_val)
	d_samples = problem.sample_d(N_val)
	y_samples = [problem.eta(theta_samples[j], d_samples[j]) for j in range(N_val)]
	
	###Use the savedata to train a series of GMMs, save them if they don't already exist; or, load them if they do
	list_gmm = []
	for N in N_list:
		gmmfile_n = "BN_model_" + str(N) + '.pkl'
		if os.path.exists(gmmfile_n):
			###Load it
			if doPrint:
				print("Loading",gmmfile_n,"...",flush=True)
			gmm_n = bn_load_gmm(gmmfile_n)
			list_gmm.append(gmm_n)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",gmmfile_n,"...",flush=True)
			gmm_n = bn_train_from_file(problem, savefile=big_savefile, do_subset=N, doPrint=True)
			
			bn_save_gmm(gmm_n, gmm_file=gmmfile_n)
			list_gmm.append(gmm_n)
	
	###Now, for each GMM, evaluate Q for each datum in the validation set
	if doPrint:
		print("Evaluating validation samples...",flush=True)
	v = [y+d for y,d in zip(y_samples, d_samples)] #turn into array
	v_traces = [None] * N_val #list of length N_val, each element is a list of length len(N_list)
	for i,vi in enumerate(v):
		if doPrint:
			print("Model evaluation",i+1,'/',len(v),"\t\t",end='\r',flush=True)
		trace = []
		for n_gmm in list_gmm:
			Qvar = gbi_var_of_conditional_pp(n_gmm, vi, verbose=0)#turn into array
			trace.append(Qvar)
		v_traces[i] = trace
	
	###Plot them all, shifting the largest-N estimate of Q to zero to hopefully see many lines converging
	v_traces_adjusted = [[abs(Qvar_n - trace[-1]) for Qvar_n in trace] for trace in v_traces]
	plt.xlabel("N for training GMM")
	plt.ylabel("Difference from final Var[Q] estimate")
	plt.xscale('log')
	for trace in v_traces_adjusted:
		plt.plot(N_list, trace, c='gray', linestyle='dashed')
	plt.show()
	
def bn_plot_data_density(problem, datafile, gmmfile, do_subset=0, doGMMPlot=False, doPrint=True):
	names = ["Q"] + problem.y_names + problem.d_names
	
	###Load data file
	qoi_train, y_d_train = bn_load_samples(problem, datafile, doPrint, do_subset=do_subset)
	training_data = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
	
	###Load gmm file
	if doGMMPlot:
		gmm = bn_load_gmm(gmmfile)
	
	###Iterate over the qyd elements, plot all of them 1d:
	for c,name in enumerate(names):
		column_data = [qyd[c] for qyd in training_data]
		print(name)
		uncertainty_prop_plot(column_data, xlab=name, c='#21aad3', saveFig='bn_plot_data_density'+name, rescaled=False, vline=[])
	
def bn_compare_model_covariance(problem, datafile, gmmfile, doPrint=True):
	###Load data file
	qoi_train, y_d_train = bn_load_samples(problem, datafile, doPrint)
	#y_d_train = [[yd[0], yd[1]] for yd in y_d_train] #stupid cut for speed
	
	###look at the covariance matrix
	data_together = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
	names_together = ["Q"] + problem.y_names + problem.d_names
	cov = covmatrix_heatmap(data_together, names_together, rescale=True) #gotta rescale to match gmm covariance!
	
	###Load gmm file
	gmm = bn_load_gmm(gmmfile)
	
	###Calculate and plot model covariance
	#https://math.stackexchange.com/questions/195911/calculation-of-the-covariance-of-gaussian-mixtures
	conditional_var_sum = sum([a*C for C,a in zip(gmm.covariances_,gmm.weights_)])
	mu_bar = sum([a * mu for mu,a in zip(gmm.means_,gmm.weights_)])
	var_conditional_mean_sum = sum([a * np.array(mu - mu_bar) * np.array(mu - mu_bar).T for mu,a in zip(gmm.means_,gmm.weights_)])
	total_cov = conditional_var_sum + var_conditional_mean_sum
	
	cov_heatmap(total_cov, names_together)
	
	###Conduct some numerical comparison
	#Frobenius norm: Measures the overall difference between matrices by summing the squared elements of the difference matrix, providing a good general measure of convergence.
	#Operator norm: Measures the maximum linear transformation induced by the difference matrix, useful when interested in the largest potential error in a specific direction.
	#Spectral norm: Measures the largest eigenvalue of the difference matrix, important when considering the impact on principal component analysis.
	#punt it for now
	0
	
def bn_evaluate_model_likelihood(problem, gmmfile, datafile="", samples=[], do_subset=0, returnSet=False, doPrint=True):
	if doPrint:
		print("Evaluating model likelihood for",gmmfile,"GMM...",flush=True)
	
	###Get validation set
	if datafile!="" and not samples:
		#use the provided training samples for evaluation
		#this is for determining convergence of the model
		if doPrint:
			print("Loading",datafile,"for training samples...",flush=True)
		qoi_train, y_d_train = bn_load_samples(problem, datafile, doPrint=True, do_subset=do_subset)
		samples = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
	elif samples and do_subset==0:
		0
	else:
		print("Insufficient instructions for bn_evaluate_model_likelihood")
		sys.exit()
		
	###Load the model to be evaluated
	gmm = bn_load_gmm(gmmfile)
		
	###Standardize the validation set
	gmm_means = gmm.standardized_mean
	gmm_stds = gmm.standardized_std
	standardized_samples = [[(yj-gmm_means[j])/gmm_stds[j] for j,yj in enumerate(y)] for y in samples]
	
	###Calculate the log likelihood of the data
	sample_loglikelihoods = gmm.score_samples(standardized_samples)
	avg_loglikelihood = np.mean(sample_loglikelihoods)
	
	if returnSet:
		return sample_loglikelihoods

	###print and return
	if doPrint:
		context_str = "validation samples" if datafile=="" else "training set samples"
		print("BN model",gmmfile,"has an average log-likelihood of",avg_loglikelihood,"for",len(standardized_samples),context_str,flush=True)
	return avg_loglikelihood

def bn_measure_likelihood_convergence(problem, big_savefile, doPrint=True):	
	N_list = [1000,4000,7000,10000,40000,70000,100000,400000,700000,1000000]
	
	###Use the savedata to train a series of GMMs, save them if they don't already exist; or, load them if they do
	list_gmm = []
	for N in N_list:
		gmmfile_n = "BN_model_" + str(N) + '.pkl'
		if os.path.exists(gmmfile_n):
			###Add its filename to list
			if doPrint:
				print("Found",gmmfile_n,"...",flush=True)
			list_gmm.append(gmmfile_n)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",gmmfile_n,"...",flush=True)
			gmm_n = bn_train_from_file(problem, savefile=big_savefile, do_subset=N, doPrint=True)
			
			bn_save_gmm(gmm_n, gmm_file=gmmfile_n)
			list_gmm.append(gmmfile_n)
	
	scores = [None] * len(N_list)
	###Now, for each GMM, evaluate the average likelihood of *its* training data
	if doPrint:
		print("Evaluating training data likelihood...",flush=True)
	for i,N in enumerate(N_list):
		#gmmfile: use BN_model_N.pkl
		#datafile: use the big huge one
		#N_val: we're not evaluating validation samples
		#do_subset: N, so that we're evaluating the N gmm with the N values used to train it
		score = bn_evaluate_model_likelihood(problem, gmmfile=list_gmm[i], datafile=big_savefile, do_subset=N, doPrint=True)
		scores[i] = val_score
	
	###Plot them all, shifting the largest-N estimate of Q to zero to hopefully see many lines converging
	plt.xlabel("N for training GMM")
	plt.ylabel("Average log-likelihood of the model's training data")
	plt.xscale('log')
	plt.plot(N_list, scores, c='r')
	plt.show()
	
def bn_measure_validation_convergence(problem, big_savefile, ncomp=0, N_list=[], N_val=0, doPrint=True, doPlot=True):	
	N_val = N_val if N_val>0 else 5		
	N_list = [1000,4000,10000,40000,100000,400000,1000000,1635000] if not N_list else N_list
	
	###Get validation set
	#this is for evaluating accuracy of the model
	if doPrint:
		print("Drawing",N_val,"validation samples...",flush=True)
	theta_val = problem.prior_rvs(N_val)
	d_val = problem.sample_d(N_val)
	y_val = [problem.eta(theta, d) for theta,d in zip(theta_val,d_val)]
	q_val = [problem.H(theta) for theta in theta_val]
	samples = [[q_val[i]]+y_val[i]+d_val[i] for i in range(N_val)]
	
	###Use the savedata to train a series of GMMs, save them if they don't already exist; or, load them if they do
	list_gmm = []
	for N in N_list:
		gmmfile_n = "BN_model_" + str(N) + '.pkl'
		if os.path.exists(gmmfile_n):
			###Add its filename to list
			if doPrint:
				print("Found",gmmfile_n,"...",flush=True)
			list_gmm.append(gmmfile_n)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",gmmfile_n,"...",flush=True)
			gmm_n = bn_train_from_file(problem, savefile=big_savefile, ncomp=ncomp, do_subset=N, doPrint=True)
			
			bn_save_gmm(gmm_n, gmm_file=gmmfile_n)
			list_gmm.append(gmmfile_n)
	
	scores = [None] * len(N_list)
	###Now, for each GMM, evaluate the average likelihood of *its* training data
	if doPrint:
		print("Evaluating training data likelihood...",flush=True)
	for i,N in enumerate(N_list):
		#gmmfile: use BN_model_N.pkl
		#datafile: we're not evaluating training samples
		#samples: validation set
		#do_subset: we're not evaluating training samples, so its irrelevant
		score = bn_evaluate_model_likelihood(problem, gmmfile=list_gmm[i], datafile=big_savefile, returnSet=False, do_subset=N, doPrint=True)
		scores[i] = score

	###Plot them all, shifting the largest-N estimate of Q to zero to hopefully see a line converging
	if doPlot:
		plt.xlabel("N for training GMM")
		plt.ylabel("Average log-likelihood of the validation set (Nv="+str()+')')
		plt.xscale('log')
		#N_plot = N_val if N_val < 100 else 100
		#for i in range(N_plot):
		#	plt.plot(N_list, [v_score[i] for v_score in scores], c='r')
		plt.plot(N_list, scores, c='r')
		plt.show()
	
#Similar to bn_train_from_file
#Except I want to evaluate BIC and MSE, and MAE for each ncomp
def bn_train_evaluate_ncomp(problem, trainfile, doPrint=True, doPlot=True):
	###Setup
	do_subset=0
	ncomps = [10,20,30,40,50,60,70,80,90,100,110,130,200]
	print("Evaluating",trainfile,"training data for GMM with number of components:",ncomps,flush=True)
	BICs = [None]*len(ncomps)
	scores = [None]*len(ncomps)

	###Load training data file
	qoi_train, y_d_train = bn_load_samples(problem, trainfile, doPrint, do_subset)
	p_mean = np.mean(qoi_train)
	p_std = np.std(qoi_train)
	yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_train]
	d_mean = np.mean(y_d_train, axis=0)
	d_std = np.std(y_d_train, axis=0)
	yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_d_train)]
	data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_train)])
	
	###Load validation data file
	#qoi_val, y_d_val = bn_load_samples(problem, valfile, doPrint, do_subset)
	#val_samples = [[qoi_val[i]]+y_d_val[i] for i in range(len(qoi_val))]
	
	###For each ncomp:
	for i,ncomp in enumerate(ncomps):
		###Train and save a model with the full data with that ncomp
		modelsave = "BN_model_"+str(len(qoi_train))+'_ncomp'+str(ncomp) +'.pkl'
		if os.path.exists(modelsave):
			###Load it
			if doPrint:
				print("Loading",modelsave,"...",flush=True)
			gmm_ncomp = bn_load_gmm(modelsave)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",modelsave,"...",flush=True)
			gmm_ncomp = gbi_train_model(qoi_train, y_d_train, verbose=2, ncomp=ncomp, careful=True)
			bn_save_gmm(gmm_ncomp, gmm_file=modelsave)
		
		###save the BICs
		BIC = gmm_ncomp.bic(data)
		BICs[i] = BIC
		print(BIC,flush=True)

	if doPlot:
		bn_train_evaluate_ncomp_plot(ncomps, BICs)
		
def bn_train_evaluate_ncomp_plot(ncomps, BICs):
	#BICs=[]
	#LRTs=[]
	print(BICs,flush=True)
	print(BICs,flush=True)
	###Extend with saved data
	#BICs.extend([])
	#ncomps = [1,10,20,30,40,50,60,70,80,90,100,110,130,200]
	
	###plot
	fig, ax1 = plt.subplots()

	# Plot the first data set on the left y-axis
	ax1.plot(ncomps, BICs, color='blue')
	ax1.set_xlabel("Number of GMM components")
	ax1.set_xticks(ncomps)
	ax1.set_ylabel('BIC', color='blue')
	ax1.tick_params(axis='y', labelcolor='blue')

	# Create the second axes sharing the x-axis
	#ax2 = ax1.twinx()

	# Plot the second data set on the right y-axis
	#ax2.plot(ncomps, LRTs, color='orange')
	#ax2.set_ylabel('LRT', color='orange')
	#ax2.tick_params(axis='y', labelcolor='orange')
	plt.show()

def bn_train_evaluate_ncomp_sanitycheck(problem, trainfile, valfile, use_val=False, doPrint=True, doPlot=True):
	###Setup
	do_subset=0
	ncomps = [1,10,20,30,40,50,60,70,80,90,100,110,130]
	print("Evaluating",trainfile,"training data for GMM with number of components:",ncomps,flush=True)
	BICs = [None]*len(ncomps)
	scores = [None]*len(ncomps)
	penalties = [None]*len(ncomps)
	
	###Load training data file
	qoi_train, y_d_train = bn_load_samples(problem, trainfile, doPrint, do_subset)
	
	if use_val:
		###Load validation data file
		qoi_data, y_d_data = bn_load_samples(problem, valfile, doPrint, do_subset)		
	else:
		###use training data
		qoi_data = qoi_train
		y_d_data = y_d_train
	
	###For each ncomp:
	for i,ncomp in enumerate(ncomps):
		###Train and save a model with the full data with that ncomp
		modelsave = "BN_model_"+str(len(qoi_train))+'_ncomp'+str(ncomp) +'.pkl'
		if os.path.exists(modelsave):
			###Load it
			if doPrint:
				print("Loading",modelsave,"...",flush=True)
			gmm_ncomp = bn_load_gmm(modelsave)
		else:
			###Train and save it
			print("Model",modelsave,"doesn't exist.",flush=True)
			sys.exit()
		
		###save the BICs - VAL
		p_mean = np.mean(qoi_data)
		p_std = np.std(qoi_data)
		yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_data]
		d_mean = np.mean(y_d_data, axis=0)
		d_std = np.std(y_d_data, axis=0)
		yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_d_data)]
		data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_data)])
		
		BIC = gmm_ncomp.bic(data)
		BICs[i] = BIC
		score = -2 * gmm_ncomp.score(data) * data.shape[0]
		scores[i] = score
		penalty = gmm_ncomp._n_parameters() * np.log(data.shape[0])
		penalties[i] = penalty
		print(BIC,gmm_ncomp.score(data),gmm_ncomp._n_parameters(),flush=True)
		#_n_parameters calculates K*((D*D - D)/2 + 2D + 1)-1, for D=62

	if doPlot:
		fig, ax1 = plt.subplots()

		# Plot the first data set on the left y-axis
		ax1.plot(ncomps, BICs, color='blue')
		ax1.set_xlabel("Number of GMM components")
		ax1.set_xticks(ncomps)
		ax1.set_ylabel('BIC', color='blue')
		ax1.tick_params(axis='y', labelcolor='blue')
		ax1.plot(ncomps, penalties, color='red')
		ax1.plot(ncomps, scores, color='green')
		plt.show()




#Do bootstrap sampling & evaluation to find 95% confidence intervals
def bn_train_convergence_confidence(problem, trainfile, N_bootstrap, ncomp, startnum=0, doPrint=True):
	do_subset=0
	###Load training data file
	qoi_train, y_d_train = bn_load_samples(problem, trainfile, doPrint, do_subset)
	
	###Iterate
	i_list = range(startnum, N_bootstrap)
	for i in i_list:
		###Make bootstrap sample
		i_samples = np.random.randint(0, len(qoi_train), size=len(qoi_train))
		qoi_bootstrap = [qoi_train[j] for j in i_samples]
		y_d_bootstrap = [y_d_train[j] for j in i_samples]
		
		###Train and save GMM on that
		modelsave = "BN_model_"+str(len(qoi_train))+'_bootstrap'+str(i)+'.pkl'
		if os.path.exists(modelsave):
			print("GMM",modelsave,"already exists, terminating.")
			sys.exit()
		else:
			if doPrint:
				print("Training and saving",modelsave,"...",flush=True)
			gmm_i = gbi_train_model(qoi_bootstrap, y_d_bootstrap, verbose=2, ncomp=ncomp, careful=True)
			bn_save_gmm(gmm_i, gmm_file=modelsave)
		
def bn_evaluate_convergence_confidence(problem, valfile, N_bootstrap, doPrint=True):
	###Load validation data file
	qoi_val, y_d_val = bn_load_samples(problem, valfile, doPrint, do_subset)

	p_mean = np.mean(qoi_bootstrap)
	p_std = np.std(qoi_bootstrap)
	yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_bootstrap]
	d_mean = np.mean(y_d_bootstrap, axis=0)
	d_std = np.std(y_d_bootstrap, axis=0)
	yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_d_bootstrap)]
	val_data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_bootstrap)])

	###Iterate
	bootstrapped_scores = [None]*N_bootstrap
	for i in range(N_bootstrap):
		###Load GMM
		modelsave = "BN_model_"+str(len(qoi_train))+'_bootstrap'+str(i)+'.pkl'
		if os.path.exists(modelsave):
			###Load it
			if doPrint:
				print("Loading",modelsave,"...",flush=True)
			gmm_g1 = bn_load_gmm(modelsave)
		else:
			print("error")
			sys.exit()
		
		###Calculate score (average likelihood of validation set), save
		gmm_i.score(val_data)
	
	##Do statistics on the scores
	#Mean
	score_mean = np.mean(bootstrapped_scores)
	
	#stddev
	score_stddev = np.std(bootstrapped_scores)
	
	#conf interval
	def conf_interval(data, conf_level):
		#returns a tuple of the low and high bounds
		return scipy.stats.t.interval(confidence=conf_level, df=len(data)-1, loc=np.mean(data), scale=scipy.stats.sem(data))
		
	conf_intervals = conf_interval(bootstrapped_scores,0.95)
	conf_interval_width = conf_intervals[1] - conf_intervals[0]
	
	return score_mean, score_stddev, conf_interval_width

