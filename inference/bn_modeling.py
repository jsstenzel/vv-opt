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
def bn_load_samples(problem, savefile, doPrint=False, doDiagnostic=False):
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
	
	#TODO do some checking in here to make sure that each line is valid, and skip if not
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

	return Q, yd

#Call bn_load_samples and gbi_train_model
def bn_train_from_file(problem, savefile, doPrint=False):
	#Load file
	qoi_train, y_d_train = bn_load_samples(problem, savefile, doPrint)
	
	#Train model	
	gmm = gbi_train_model(qoi_train, y_d_train, verbose=2, ncomp=0)
	
	#Print and return
	if doPrint:
		print("Trained GMM with",len(qoi_train),"samples from",savefile)
		print(gmm, flush=True)
		
	return gmm
	
#This calculates MSE between the truth and model
#Note that of course, conditioning the GMM gets you the posterior predictive of Q
#So what I actually want to compare here is how well Var[GMM(y,d)] compares to Var[Q]
#NOTE that this doesn't compare Var[GMM(y,d)] to Var[p(Q|y)], to do that I need to do inference, such as with MCMC
def bn_measure_model_mse(problem, gmm, N, doPrint=True):
	###Draw N totally new validation samples
	theta_samples = problem.prior_rvs(N)
	d_samples = problem.sample_d(N)
	
	###Evaluate the validation set for "truth" case
	#for each theta,d get a bunch of subsamples of Q, and calculate the variance of that
	subsample_N = N #i suppose?
	Qvar_true = []
	for i in range(N):
		qoi_train = [problem.H(theta_samples[i]) for i in range(subsample_N)]
		#y_train = [problem.eta(theta_samples[i], d_samples[i]) for _ in range(subsample_N)]
		Qvar_true.append(np.var(qoi_train, ddof=0))  #not ddof=1, since I want to compare with the conditioned GMM, which uses the "actual" variance!
	
	###Evaluate the validation set with the GMM
	y_samples = [problem.eta(theta_samples[j], d_samples[j]) for j in range(N)]
	#this introduces uncertainty that is not present in the truth case. To calculate real MSE, i would need to actually calculate posterior
		
	v = [y+d for y,d in zip(y_samples, d_samples)] #turn into array
	Qvar_gmm = [gbi_var_of_conditional_pp(gmm, vi) for vi in v] #turn into array
	
	SE = [(truth - model)**2 for truth,model in zip(Qvar_true, Qvar_model)]
	MSE = np.mean(SE)
	if doPrint:
		print("MSE for",N,"validation samples:",MSE, flush=True)
		
	return MSE

#Save the GMM to file for easy grabbing later
def bn_save_gmm(gmm, gmm_file):
	ncomp = gmm.n_components
	weights = gmm.weights_
	means = gmm.means_
	covariances = gmm.covariances_
	
	filename = gmm_file if gmm_file.endswith('.csv') else gmm_file+'.csv'
	with open(filename, 'w+', newline='') as csvfile:
		writer = csv.writer(csvfile)
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
			weight_lists.append(float(row[0]))
			mean_lists.append(float(row[1]))
			cov_lists.append([float(e) for e in row[2:]])
	
	#make the parameter arrays
	weights = np.array(weight_lists)
	ncomp = len(weights)
	
	means = np.array(mean_lists) #ncomp x d
	d = len(means[0])
	covariances = [np.array(cov).reshape(d,d) for cov in cov_lists] #ncomp x dxd #reshape each list of floats into a dxd array
	
	#ugly but i gotta
	#https://stackoverflow.com/questions/68541971/sklearn-sample-from-gaussianmixture-without-fitting
	data_gmm = GaussianMixture(n_components=ncomp)
	data_gmm.weights_ = weights
	data_gmm.means_ = means
	data_gmm.covariances_ = covariances

	return data_gmm