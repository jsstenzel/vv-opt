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
from sklearn.mixture import GaussianMixture

sys.path.append('..')
from inference.goal_based_inference import *
from uq.plotmatrix import *

#Draw N samples and save in a consistent 1-line format to savefile
#Do it iteratively, to support parallelization and clustering
def bn_sampling(problem, savefile, N, buffer_rate=1, doPrint=False):
	filename = savefile if savefile.endswith('.csv') else savefile+'.csv'
	
	data_buffer = []
	for i in range(N):
		if doPrint:
			print("Drawing Q & y samples",i,"...",flush=True)
		#Drawing samples theta, d
		theta_sample = problem.prior_rvs(1)
		d_sample = problem.sample_d(1)
			
		#Model propagation for Q, y
		qoi_train = problem.H(theta_sample)
		y_train = problem.eta(theta_sample, d_sample)
		
		#Append the new BN sample to file
		save_data = y_train + d_sample + [qoi_train]
		
		if buffer_rate <= 1: #just stream it into the save file
			with open(filename, 'a+', newline='') as csvfile:
				writer = csv.writer(csvfile)
				writer.writerow(save_data)
		else: #save to buffer, and dump buffer to file at buffer_rate
			data_buffer.append(save_data)
			if (i+1) % buffer_rate == 0:
				print("(dump)",flush=True)
				with open(filename, 'a+', newline='') as csvfile:
					writer = csv.writer(csvfile)
					for data in data_buffer:
						writer.writerow(data)
				data_buffer.clear()
				

#Read and interpret the savefile
def bn_load_samples(problem, savefile, doPrint=False, do_subset=0, doDiagnostic=False):
	filename = savefile if savefile.endswith('.csv') else savefile+'.csv'

	###Make sure the files exist
	if not os.path.isfile(filename):
		print("File",filename,"is missing")
		sys.exit()
		
	###Safely read out all of the samples into matrices
	y = []
	d = []
	Q = []
	
	if doPrint:
		print("Reading the data files...",flush=True)
	
	with open(filename) as csvfile:
		csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
		for l,row in enumerate(csvreader):
			if len(row) != problem.dim_y + problem.dim_d + 1:
				if doDiagnostic:
					print("Warning: dropped line",l+1,"(length "+str(len(row))+' expected', str(problem.dim_y + problem.dim_d + 1)+')',"from",filename)
			else:
				ygrab = [float(e) for e in row[:problem.dim_y]]
				dgrab = [float(e) for e in row[problem.dim_y:-1]]
				Qgrab = float(row[-1])
				y.append(ygrab) #should be length dim_y
				d.append(dgrab) #should be length dim_d = row - dim_y - 1
				Q.append(Qgrab)
	
	#zip it together at the end, in case i ever need them separate for something later
	yd = [y_i + d_i for y_i,d_i in zip(y,d)]
	
	if do_subset:
		Q = Q[:do_subset]
		yd = yd[:do_subset]

	return Q, yd

#Call bn_load_samples and gbi_train_model
def bn_train_from_file(problem, savefile, do_subset=0, doPrint=False):
	###Load file
	qoi_train, y_d_train = bn_load_samples(problem, savefile, doPrint, do_subset)
	#y_d_train = [[yd[0], yd[1]] for yd in y_d_train] #stupid cut for speed
	
	###look at the covariance matrix and eigenvalues of the data, see if we have bad correlation
	if False:
		data_together = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
		names_together = ["Q"] + problem.y_names + problem.d_names
		cov = covmatrix_heatmap(data_together, names_together, rescale=False)
		eigenval, _ = np.linalg.eig(cov)
		print("eigenvalues of the Qyd covariance matrix:",eigenval)
	
	###Train model	
	gmm = gbi_train_model(qoi_train, y_d_train, verbose=2, ncomp=0, careful=True)
	
	###Print and return
	if doPrint:
		print("Trained GMM with",len(qoi_train),"samples from",savefile)
		print(gmm,len(gmm.means_[0]),"dimensions", flush=True)
		
	return gmm
	
#This calculates MSE between the truth and model
#Note that of course, conditioning the GMM gets you the posterior predictive of Q
#So what I actually want to compare here is how well Var[GMM(y,d)] compares to Var[Q]
#NOTE that this doesn't compare Var[GMM(y,d)] to Var[p(Q|y)], to do that I need to do inference, such as with MCMC
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

#Save the GMM to file for easy grabbing later
def bn_save_gmm(gmm, gmm_file):
	ncomp = gmm.n_components
	weights = gmm.weights_
	means = gmm.means_
	covariances = gmm.covariances_
	norm_means = gmm.standardized_mean
	norm_stds = gmm.standardized_std
	
	filename = gmm_file if gmm_file.endswith('.csv') else gmm_file+'.csv'
	with open(filename, 'w+', newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(norm_means)
		writer.writerow(norm_stds)
		for i in range(ncomp):
			covlist = covariances[i].flatten().tolist()
			newrow = [weights[i]] + means[i].tolist() + covlist
			writer.writerow(newrow)
			
	print("GMM saved to",filename,flush=True)

#Load GMM from file
def bn_load_gmm(gmm_file):
	weight_lists = []
	mean_lists = []
	cov_lists = []

	filename = gmm_file if gmm_file.endswith('.csv') else gmm_file+'.csv'
	with open(filename) as csvfile:
		csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
		#each row represents one of the Gaussian components
		for l,row in enumerate(csvreader):
			if l==0:
				norm_means = [float(mu_i) for mu_i in row]
			elif l==1:
				norm_stds = [float(s_i) for s_i in row]
			else:
				param_len = int(max(np.roots([1,1,1-len(row)]))) #lol line format is 1 + x + x^2, find x=# params
				weight_lists.append(float(row[0]))
				mean_lists.append([float(x) for x in row[1:param_len+1]])
				cov_lists.append([float(e) for e in row[param_len+1:]])
	
	#make the parameter arrays
	weights = np.array(weight_lists)
	ncomp = len(weights)
	
	means = np.array(mean_lists) #ncomp x d
	d = len(means[0])
	try:
		covariances = [np.array(cov).reshape(d,d) for cov in cov_lists] #ncomp x dxd #reshape each list of floats into a dxd array
	except:
		print("Error:",gmm_file,"does not have the right format of weight, mean, and covariance on each row")
		print("file row lengths:")
		for i in range(ncomp):
			print(1 + len(mean_lists[i]) + len(cov_lists[i]))
		sys.exit()
	
	#ugly but i gotta
	#https://stackoverflow.com/questions/68541971/sklearn-sample-from-gaussianmixture-without-fitting
	data_gmm = GaussianMixture(n_components=ncomp)
	data_gmm.weights_ = weights
	data_gmm.means_ = means
	data_gmm.covariances_ = covariances
	data_gmm.standardized_mean = norm_means
	data_gmm.standardized_std = norm_stds

	return data_gmm
	

def bn_measure_stability_convergence(problem, big_savefile, N_val, doPrint=True):
	N_val = N_val if N_val>0 else 5		
	N_list = [1000,4000,10000,40000,100000,400000]
	
	###Draw N_val totally new validation samples
	if doPrint:
		print("Drawing",N_val,"validation samples...",flush=True)
	theta_samples = problem.prior_rvs(N_val)
	d_samples = problem.sample_d(N_val)
	y_samples = [problem.eta(theta_samples[j], d_samples[j]) for j in range(N_val)]
	
	###Use the savedata to train a series of GMMs, save them if they don't already exist; or, load them if they do
	list_gmm = []
	for N in N_list:
		gmmfile_n = "BN_model_" + str(N) + '.csv'
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
	
def bn_evaluate_model_likelihood(problem, gmmfile, datafile="", N_val=0, do_subset=0, doPrint=True):
	if doPrint:
		print("Evaluating model likelihood for",gmmfile,"GMM...",flush=True)
	
	###Get validation set
	if datafile=="" and N_val>0:
		#calculate some validation samples
		#this is for evaluating accuracy of the model
		if doPrint:
			print("Drawing",N_val,"validation samples...",flush=True)
		theta_val = problem.prior_rvs(N_val)
		d_val = problem.sample_d(N_val)
		y_val = [problem.eta(theta, d) for theta,d in zip(theta_val,d_val)]
		q_val = [problem.H(theta) for theta in theta_val]
		samples = [[q_val[i]]+y_val[i]+d_val[i] for i in range(N_val)]
	elif datafile!="":
		#use the provided training samples for evaluation
		#this is for determining convergence of the model
		if doPrint:
			print("Loading",datafile,"for training samples...",flush=True)
		qoi_train, y_d_train = bn_load_samples(problem, datafile, doPrint=True, do_subset=do_subset)
		samples = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
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

	###print and return
	if doPrint:
		context_str = "validation samples" if datafile=="" else "training set samples"
		print("BN model",gmmfile,"has an average log-likelihood of",avg_loglikelihood,"for",len(standardized_samples),context_str,flush=True)
	return avg_loglikelihood


#def bn_measure_likelihood_convergence(problem, datafile="", gmmfile, doPrint=True):