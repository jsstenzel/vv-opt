#import argparse
import os
import sys
import csv

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from itertools import islice

sys.path.insert(0, "..")
from uq.uncertainty_propagation import *
from uq.plotmatrix import *
#from problems.problem_definition import *

#Here, I am manually implementing Sobol Saltelli algorithm

"""
Summary of functions:

problem_saltelli_sample
	Runs the full sample + eval + indices problem for N samples
	Will use saved samples or create more as necessary
	Calls saltelli_eval_sample, then saltelli_indices
	
problem_saltelli_sobol
	Runs the full sample + eval + indices problem for N samples
	Will use saved samples or create more as necessary
	Calls saltelli_eval_sobol, then saltelli_indices

saltelli_eval
	Take in already-made samples
	construct new A, B, Ci
	evaluate model for all new samples
	append new samples & evaluations onto files
	
saltelli_indices
	Take in filenames
	Calculate and print Si, STi
		
saltelli_eval_sample
	Take in variable names, distributions, bounds/parameters
	Samples from those distributions
	Then calls saltelli_eval

saltelli_eval_sobol
	Take in variable names, distributions, bounds/parameters
	Samples from the SciPy Sobol sequence, starting the sequence in a place that continues any existing Sobol samples
	Performs CDF post-processing on the Sobol sequence to make it match the correct distributions
	Then calls saltelli_eval
	
Example: Saltelli iteration over M
	Take in a list of cumulative sample sizes to use
	Take in a sampler function with 1 argument: sample size
	Function loop:
		Go to each provided sample size
		use saltelli_eval to make up the samples, if necessary
		use saltelli_indices
"""

def problem_saltelli_sample(N, base_name, var_names, var_dists, var_params, model, doPrint=True):
	###First, figure out how many new samples are actually needed to get to N
	old_samples = 0
	if not os.path.isfile(base_name+'_A.csv') or not os.path.isfile(base_name+'_B.csv') or any([not os.path.isfile(base_name+'_C_'+name+'.csv') for name in var_names]):
		#clear_saltelli_sample_files(base_name, var_names) #in case of weird partial files
		0
	else:
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for row in csvreader:
				if any(row):  # Check if any of the fields in the row are non-empty
					old_samples += 1
	old_samples *= 2 #because there are different samples located in both _A and _B
		
	if doPrint:
		print("problem_saltelli_sample for N=", N, "running on",base_name,"which has",old_samples,"old samples.",flush=True)

	if old_samples < N:
		###Get the new samples and return indices
		new_samples = N - old_samples
		saltelli_eval_sample(base_name, new_samples, var_names, var_dists, var_params, model, doPrint=doPrint)
		return saltelli_indices(base_name, var_names, doPrint=doPrint)
	elif old_samples == N:
		###Return indices
		return saltelli_indices(base_name, var_names, doPrint=doPrint)
	else: #old_samples > N
		###Return indices with a subset of the samples
		return saltelli_indices(base_name, var_names, do_subset=N, doPrint=doPrint)
	
	
def problem_saltelli_sobol(p, base_name, var_names, var_dists, var_params, model, doPrint=True):
	###First, figure out how many new samples are actually needed to get to 2**p
	old_samples = 0
	if not os.path.isfile(base_name+'_A.csv') or not os.path.isfile(base_name+'_B.csv') or any([not os.path.isfile(base_name+'_C_'+name+'.csv') for name in var_names]):
		#clear_saltelli_sample_files(base_name, var_names) #in case of weird partial files
		0
	else:
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for row in csvreader:
				if any(row):  # Check if any of the fields in the row are non-empty
					old_samples += 1
	old_samples *= 2 #because there are different samples located in both _A and _B
		
	if doPrint:
		print("problem_saltelli_sobol for N=", 2**p, "running on",base_name,"which has",old_samples,"old samples.",flush=True)

	if old_samples < 2**p:
		###Get the new samples and return indices
		saltelli_eval_sobol(base_name, p, var_names, var_dists, var_params, model, doPrint=doPrint) #p is total power, not new samples
		return saltelli_indices(base_name, var_names, doPrint=doPrint)
	elif old_samples == 2**p:
		###Return indices
		return saltelli_indices(base_name, var_names, doPrint=doPrint)
	else:
		###Return indices with a subset of the samples
		return saltelli_indices(base_name, var_names, do_subset=2**p, doPrint=doPrint)


def saltelli_eval(new_samples, base_name, var_names, model, doPrint=False):
	if len(new_samples) % 2 != 0:
		print("Need an even number of new samples to saltelli_eval", flush=True)
		sys.exit()
		
	new_A = []
	new_B = []
	###Split 2M into A and B "deterministically"
	for m,row in enumerate(new_samples):
		if m % 2 == 0:
			new_A.append(list(row))
		else:
			new_B.append(list(row))
	#new_A = list(new_samples[:int(len(new_samples)/2)]) #doesnt work if its a deep np array. but i dont want to deal with checking for whether it is or not
	#new_B = list(new_samples[int(len(new_samples)/2):])
			
	###Construct the p Ci matrices
	new_C = []
	for p,name in enumerate(var_names):
		new_Ci = []
		
		#Construct new_Ci by subsituting the pth column of A into B
		for n,b_row in enumerate(new_B):
			a_row = new_A[n]
			
			c_row = b_row[:]
			c_row[p] = a_row[p]
			new_Ci.append(c_row)
		
		#Now save that into a p-long list of C's
		new_C.append(new_Ci)
		
	###Perform the M*(2+p) model evaluations, append onto each row of the new matrices
	if doPrint:
		nevals = int((len(new_samples)/2) * (2+len(var_names)))
		print("Performing",nevals,"model evaluations ...",flush=True)
		
	for m,a_row in enumerate(new_A):
		yA = model(a_row)
		new_A[m].append(yA)
	
	for m,b_row in enumerate(new_B):
		yB = model(b_row)
		new_B[m].append(yB)
	
	for new_Ci in new_C:
		for m,ci_row in enumerate(new_Ci):
			yCi = model(ci_row)
			new_Ci[m].append(yCi)
		
	###Finally, append these new matrices onto each of the 2+p files
	if doPrint:
		print("Saving samples to",base_name,"...",flush=True)
		
	aw = 'a+'# if os.path.exists(base_name+'_A.csv') else 'w+'
	with open(base_name+'_A.csv', aw, newline='') as csvfile:
		writer = csv.writer(csvfile)
		for row in new_A:
			writer.writerow(row)
	
	aw = 'a+'# if os.path.exists(base_name+'_B.csv') else 'w+'
	with open(base_name+'_B.csv', aw, newline='') as csvfile:
		writer = csv.writer(csvfile)
		for row in new_B:
			writer.writerow(row)

	for p,new_Ci in enumerate(new_C):
		aw = 'a+'# if os.path.exists(base_name+'_C_'+var_names[p]+'.csv') else 'w+'
		with open(base_name+'_C_'+var_names[p]+'.csv', aw, newline='') as csvfile:
			writer = csv.writer(csvfile)
			for row in new_Ci:
				writer.writerow(row)

#do_subset argument lets us use a subset of the saved samples. Should be <= 2M
def saltelli_indices(base_name, var_names, do_subset=0, doPrint=True):
	if not os.path.isfile(base_name+'_A.csv'):
		print("File",base_name+'_A.csv',"is missing")
		sys.exit()
	if not os.path.isfile(base_name+'_B.csv'):
		print("File",base_name+'_B.csv',"is missing")
		sys.exit()
	for p,name in enumerate(var_names):
		if not os.path.isfile(base_name+'_C_'+name+'.csv'):
			print("File",base_name+'_C_'+name+'.csv',"is missing")
			sys.exit()

	###First, read out all of the samples into matrices
	Ay = []
	By = []
	Cy = []
	
	if doPrint:
		print("Reading the data files...",flush=True)
	
	if do_subset == 0:
		with open(base_name+'_A.csv', errors='ignore') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for row in csvreader:
				Ay.append([float(elem) for elem in row])
		#This change assumes that the data files are well-behaved; use SA_datalist_health_check first
		M = len(Ay)
		pp = len(var_names)
		By = np.zeros((M,pp+1))
		Cy = []
		
		with open(base_name+'_B.csv', errors='ignore') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for i,row in enumerate(csvreader):
				for e,elem in enumerate(row):
					By[i][e] = float(elem)

		for p,name in enumerate(var_names):
			Ciy = np.zeros((M,pp+1))
			with open(base_name+'_C_'+name+'.csv', errors='ignore') as csvfile:
				csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
				for i,row in enumerate(csvreader):
					for e,elem in enumerate(row):
						Ciy[i][e] = float(elem)
			Cy.append(Ciy)
	else:
		lim = int(do_subset/2)
		###Optionally, we can analyze less than the full set of provided samples
		with open(base_name+'_A.csv', errors='ignore') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for row in islice(csvreader, lim):
				Ay.append([float(elem) for elem in row])
		#This change assumes that the data files are well-behaved; use SA_datalist_health_check first
		M = len(Ay)
		pp = len(var_names)
		By = np.zeros((M,pp+1))
		Cy = []
		
		with open(base_name+'_B.csv', errors='ignore') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for i,row in enumerate(islice(csvreader, lim)):
				By[i] = [float(elem) for elem in row]

		for p,name in enumerate(var_names):
			Ciy = np.zeros((M,pp+1))
			with open(base_name+'_C_'+name+'.csv', errors='ignore') as csvfile:
				csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
				for row in islice(csvreader, lim):
					for e,elem in enumerate(row):
						Ciy[i][e] = float(elem)
			Cy.append(Ciy)
		
	###Isolate A,B,Ci and yA,yB,yCi
	if doPrint:
		print("Processing samples...",flush=True)
		
	A = [Ay_row[:-1] for Ay_row in Ay] #all but last element
	yA = [Ay_row[-1] for Ay_row in Ay] #only last element
	B = [By_row[:-1] for By_row in By]
	yB = [By_row[-1] for By_row in By]
	
	C = []
	yC = []
	for p,Ciy in enumerate(Cy):
		Ci = [Ciy_row[:-1] for Ciy_row in Ciy]
		yCi = [Ciy_row[-1] for Ciy_row in Ciy]
		C.append(Ci)
		yC.append(yCi)

	if doPrint:
		print("Calculating indices...",flush=True)
	
	###Provide final calculations
	f02 = np.dot(np.mean(yA),np.mean(yA))  #inner product of p-length and p-length vectors
	S = []
	ST = []
	model_evals = M*(len(var_names)+2)
	
	yAyA = np.dot(yA, yA)      #inner product of n-length and n-length vectors
	for p,name in enumerate(var_names):
		yAyCi = np.dot(yA, yC[p])  #inner product of n-length and n-length vectors
		yByCi = np.dot(yB, yC[p])  #inner product of n-length and n-length vectors
		#print(yAyCi, yAyA, yByCi)
		Si = (yAyCi/M - f02)/(yAyA/M - f02)
		STi = 1 - (yByCi/M - f02)/(yAyA/M - f02)
		S.append(Si)
		ST.append(STi)
		
	if doPrint:
		print("-------------------------------------")
		print("Saltelli GSA algorithm, M =",M,"p =",pp)
		print("Total number of model evaluations:",model_evals)
		print("Parameter",'\t','\t',"S_i",)
		for p,name in enumerate(var_names):
			print(name,'\t','\t',S[p])
		print("-------------------------------------")
		print("Parameter",'\t','\t',"S_Ti")
		for p,name in enumerate(var_names):
			print(name,'\t','\t',ST[p])
		print("-------------------------------------", flush=True)
	
	return S, ST, model_evals
	
def saltelli_eval_sample(base_name, sample_size, var_names, var_dists, var_params, model, doPrint=True):
	###Process the distributions
	allowable_dists = {
	'gaussian': 2, #mu, sigma
	'gaussian_multivar': 2, #mean vector, covariance
	'gamma_ab': 2, #alpha, beta
	'gamma_mv': 2, #mean, variance
	'beta': 2, #a, b
	'lognorm': 2, #mean, variance
	'uniform': 2, #left, right
	'nonrandom': 1, #return value
	'gp_expquad': 4 #variance, ls, prior_pts, mean_fn
	}
	for name,dtype,params in zip(var_names,var_dists,var_params):
		if dtype not in allowable_dists.keys():
			print("Distribution ",dtype,"for",name,"not recognized by saltelli_eval_sample")
			sys.exit()
		if len(params) != allowable_dists[dtype]:
			print("Distribution ",dtype,"for",name,"expects",allowable_dists[dtype],"params, not",len(params))
			sys.exit()
			
	if doPrint:
		print("Sampling...",flush=True)
	
	###Samples from those distributions
	vals = [] #a list length sample_size of random numbers of size sample_size
	val_names = []
	for name,dtype,params in zip(var_names,var_dists,var_params): ###iterate over p
		#generate the rvs for this one particular theta
		if dtype == 'gaussian':
			mu = params[0]
			sigma = params[1]
			samples = scipy.stats.norm.rvs(size=sample_size, loc=mu, scale=sigma)
			vals.append(samples.tolist())
			val_names.append(name)
		elif dtype == 'gaussian_multivar':
			mean_vector = np.array(params[0])
			covariance = np.array(params[1])
			multisamples = scipy.stats.multivariate_normal.rvs(size=sample_size, mean=mean_vector, cov=covariance)
			i=0
			for samples in multisamples.T: #TODO check this later
				vals.append(samples.tolist())
				val_names.append(name+"_"+str(i))
				i+=1
		elif dtype == 'gamma_ab':
			alpha = params[0]
			beta = params[1]
			samples = scipy.stats.gamma.rvs(size=sample_size, a=alpha, scale=1.0/beta)
			vals.append(samples.tolist())
			val_names.append(name)
		elif dtype == 'gamma_mv':
			mean = params[0]
			variance = params[1]
			alpha = mean**2 / variance
			beta = mean / variance
			samples = scipy.stats.gamma.rvs(size=sample_size, a=alpha, scale=1.0/beta)
			vals.append(samples.tolist())
			val_names.append(name)
		elif dtype == 'beta':
			a = params[0]
			b = params[1]
			samples = scipy.stats.beta.rvs(a=a, b=b, size=sample_size)
			vals.append(samples.tolist())
			val_names.append(name)
		elif dtype == 'lognorm':
			mu = params[0]
			sigma = params[1]
			samples = scipy.stats.lognorm.rvs(size=sample_size, s=sigma, scale=np.exp(mu))
			vals.append(samples.tolist())
			val_names.append(name)
		elif dtype == 'uniform':
			left = params[0]
			right = params[1]
			samples = scipy.stats.uniform.rvs(size=sample_size, loc=left, scale=right-left) #dist is [loc, loc + scale]
			vals.append(samples.tolist())
			val_names.append(name)
		elif dtype == 'gp_expquad':
			variance = params[0]
			ls = params[1]
			prior_pts = params[2]
			mean_fn = params[3]
			gp_prior = GaussianProcessSample1D(variance, ls, prior_pts, mean_fn)
			thetas_i = [gp_prior.sample() for _ in range(num_vals)]
			vals.append(thetas_i)
			val_names.append(name)
		else:
			raise ValueError("saltelli_eval_sample did not expect prior type "+str(dtype))
			
	#turn from list of rvs at each prior, to list of theta rvs
	vals = np.transpose(vals)
	
	if doPrint:
		print("Performing evaluation...",flush=True)
	
	###Then calls saltelli_eval
	return saltelli_eval(vals, base_name, val_names, model, doPrint=doPrint)
	

#Sobol sequence only has nice properties when sample size is a power of two
#For this function to ensure that this is always the case, you don't specify the number of new points you want to sample
#Instead, you specify the total power you want the sequence to have
#Then, the necessary new samples from the Sobol sequence are drawn in a way that brings you from n_prev to 2^total_sample_power
def saltelli_eval_sobol(base_name, total_sample_power, var_names, var_dists, var_params, model, doPrint=True):
	###Process the distributions
	allowable_dists = {
	'gaussian': 2, #mu, sigma
	'gaussian_multivar': 2, #mean vector, covariance
	'gamma_ab': 2, #alpha, beta
	'gamma_mv': 2, #mean, variance
	'beta': 2, #a, b
	'lognorm': 2, #mean, variance
	'uniform': 2, #left, right
	'nonrandom': 1, #return value
	'gp_expquad': 4 #variance, ls, prior_pts, mean_fn
	}

	for name,dtype,params in zip(var_names,var_dists,var_params):
		if dtype not in allowable_dists.keys():
			print("Distribution ",dtype,"for",name,"not recognized by saltelli_eval_sample")
			sys.exit()
		elif len(params) != allowable_dists[dtype]:
			print("Distribution ",dtype,"for",name,"expects",allowable_dists[dtype],"params, not",len(params))
			sys.exit()
		elif dtype == 'gaussian_multivar':
			print("saltelli_eval_sobol can't do gaussian_multivar yet. I don't know of a good way to do inverse CDF transform for the multivariate Gaussian; scipy doesn't implement that")
			sys.exit()
		elif dtype == 'gp_expquad':
			print("saltelli_eval_sobol can't do gp_expquad, consider using hyperparameter(s) instead")
			sys.exit()
	
	###Samples from the SciPy Sobol sequence, starting the sequence in a place that continues any existing Sobol samples
	#First, find where our starting point is by looking at the already-generated sequences
	M = 0
	if os.path.isfile(base_name+'_A.csv') and os.path.isfile(base_name+'_B.csv'):
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for row in csvreader:
				if any(row):  # Check if any of the fields in the row are non-empty
					M += 1
	existing_points = 2*M #We've already used the first 2M items in the Sobol sequence, now lets start at 2M+1
	if M!=0 and existing_points & (existing_points - 1) != 0:
		print("Warning, something unexpected for saltelli_eval_sobol -- there are",existing_points,"pre-existing samples saved at",base_name,"which is not a power of 2.",flush=True)
		#This suggests you may have used the wrong files
	
	#Now, figure out how many more points to grab
	if total_sample_power <= 0: #special behavior: just go up 1 power
		p_prev = int(np.log2(existing_points))
		new_points = 2**(p_prev+1) - existing_points
	elif existing_points >= 2**total_sample_power:
		print("saltelli_eval_sobol was asked to make",2**total_sample_power,"total samples, but there are already",existing_points,"saved at",base_name)
		sys.exit()
	else:
		new_points = 2**total_sample_power - existing_points
		
	if doPrint:
		print("Sampling...",flush=True)
	
	#Now, grab the necessary Sobol sequence values	
	#TODO gotta totally rewrite this, cant really rely on random_base2
	#Gota use random(n=) in order to get from existing_points to 2^total_sample_power
	sampler = scipy.stats.qmc.Sobol(d=len(var_names), scramble=False)
	if existing_points > 0:
		sampler.fast_forward(existing_points)
	sample = sampler.random(n=new_points)
	
	###Performs CDF post-processing on the Sobol sequence to make it match the correct distributions
	#Think about this like im iterating through the columns of sample
	for j,dtype in enumerate(var_dists):
		name = var_names[j]
		params = var_params[j]
		#In each column, I find the appropriate distribution
		#Then, I run dtype.ppf(column, params), which calculates all of those ppf's independently
		#This gives me the appropriate sample value at each row, which I sub in
		if dtype == 'gaussian':
			mu = params[0]
			sigma = params[1]
			sample[:,j] = scipy.stats.norm.ppf(sample[:,j], loc=mu, scale=sigma)
		elif dtype == 'gamma_ab':
			alpha = params[0]
			beta = params[1]
			sample[:,j] = scipy.stats.gamma.ppf(sample[:,j], a=alpha, scale=1.0/beta)
		elif dtype == 'gamma_mv':
			mean = params[0]
			variance = params[1]
			alpha = mean**2 / variance
			beta = mean / variance
			sample[:,j] = scipy.stats.gamma.ppf(sample[:,j], a=alpha, scale=1.0/beta)
		elif dtype == 'beta':
			a = params[0]
			b = params[1]
			sample[:,j] = scipy.stats.beta.ppf(sample[:,j], a=a, b=b)
		elif dtype == 'lognorm':
			mu = params[0]
			sigma = params[1]
			sample[:,j] = scipy.stats.lognorm.ppf(sample[:,j], s=sigma, scale=np.exp(mu))
		elif dtype == 'uniform':
			left = params[0]
			right = params[1]
			sample[:,j] = [sample*(right-left)+left for sample in sample[:,j]]
			#TODO check on that, I don't think that's quite the right shape
		else:
			raise ValueError("saltelli_eval_sobol did not expect prior type "+str(dtype))
			
	if doPrint:
		print("Performing evaluation...",flush=True)
	
	###Then calls saltelli_eval
	return saltelli_eval(sample, base_name, var_names, model)


def clear_saltelli_sample_files(base_name, var_names):
	os.remove(base_name+'_A.csv')
	os.remove(base_name+'_B.csv')
	for p,name in enumerate(var_names):
		os.remove(base_name+'_C_'+var_names[p]+'.csv')

def __ishigami(X):
	a=7
	b=0.1
	x1 = X[0]
	x2 = X[1]
	x3 = X[2]
	return np.sin(x1) + a * np.sin(x2)**2 + b * x3**4 * np.sin(x1)
	
#According to https://uqpyproject.readthedocs.io/en/stable/auto_examples/sensitivity/comparison/ishigami.html#sphx-glr-auto-examples-sensitivity-comparison-ishigami-py
#I should expect the following first-order indices:
#S1 = 0.3139
#S2 = 0.4424
#S3 = 0.0
#and the following total-order indices:
#ST1 = 0.55758886
#ST2 = 0.44241114
#ST3 = 0.24368366

def __print_ishigami_indices(a=7, b=0.1):
	print("Expected first-order indices:")
	VY = a**2/8 + (b*np.pi**4)/5 + b**2 * (np.pi**8)/18 + 0.5 
	V1 = 0.5 * (1 + (b*np.pi**4)/5)**2
	V2 = a**2 / 8
	V3 = 0
	print("S_1 =", V1/VY)
	print("S_2 =", V2/VY)
	print("S_3 =", V3/VY)
	print("Expected total-effect indices:")
	VT1 = 0.5 * (1 + (b*np.pi**4)/5)**2 + (8*b*b*np.pi**8)/225
	VT2 = a**2 / 8
	VT3 = (8*b*b*np.pi**8)/225
	print("S_T1 =", VT1/VY)
	print("S_T2 =", VT2/VY)
	print("S_T3 =", VT3/VY)
	
def __compare_samplers(n, base1 = "ishigami", base2 = "ishigami_sobol"):
	###compare samplers
	samples1 = []
	samples2 = []
	
	with open(base1+'_A.csv') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in islice(csvreader, n):
			samples1.append([float(elem) for elem in row])

	with open(base1+'_B.csv') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in islice(csvreader, n):
			samples1.append([float(elem) for elem in row])

	with open(base2+'_A.csv') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in islice(csvreader, n):
			samples2.append([float(elem) for elem in row])

	with open(base2+'_B.csv') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in islice(csvreader, n):
			samples2.append([float(elem) for elem in row])
	
	for i in range(len(samples1[0])):
		uncertainty_prop([row[i] for row in samples1], doPlot=True, doPrint=True)
		uncertainty_prop([row[i] for row in samples2], doPlot=True, doPrint=True)
		
	covmatrix_heatmap(samples1, ["x1","x2","x3","y"])
	covmatrix_heatmap(samples2, ["x1","x2","x3","y"])
		
	sys.exit()
	
def plot_gsa_histograms(varnames, S=None, ST=None, title="", logplot=False):
	things = [S, ST]
	for i,thing in enumerate(things):
		if thing is not None:
			# set width of bar
			barWidth = 0.2
			#fig = plt.subplots(figsize =(10, 5))

			# Make the plot
			plt.bar(varnames, ST, color='orange', width = barWidth, edgecolor='grey',)
			plt.title(title)
			
			if i==0:
				plt.ylabel('First-order Sobol index', fontweight ='bold', fontsize = 10)
			elif i==1:
				plt.ylabel('Total-effect Sobol index', fontweight ='bold', fontsize = 10)
			
			# Adding Xticks
			plt.xlabel('Parameters', fontweight ='bold', fontsize = 10)
			#plt.xticks([r + barWidth for r in range(len(varnames))], varnames)
			plt.xticks(rotation=90)
			plt.tight_layout()
			if logplot:
				plt.yscale('log')	
				
			#add names (no whiskers)
			plt.errorbar(varnames, thing, yerr=0, fmt='|', color='k')
			plt.grid(axis = 'y')
			
			plt.show()

if __name__ == '__main__':  
	__compare_samplers(2**12)
	###########################
	#Full demonstration
	###########################
	
	###Set up the problem
	var_names = ['x1','x2','x3']
	dists = ['uniform', 'uniform', 'uniform']
	bounds = [[-np.pi,np.pi],[-np.pi,np.pi],[-np.pi,np.pi]]
	
	###Run two convergence tests
	list_p = [10]
	list_N = [2**p for p in list_p]
	list_S_random = []
	list_S_sobol = []
	
	for i in range(len(list_p)):
		Sr, STr, _ = problem_saltelli_sample(list_N[i], "ishigami", var_names, dists, bounds, __ishigami, doPrint=(i==len(list_p)-1))
		list_S_random.append(Sr)
		
		Ss, STs, _ = problem_saltelli_sobol(list_p[i], "ishigami_sobol", var_names, dists, bounds, __ishigami, doPrint=(i==len(list_p)-1))
		list_S_sobol.append(Ss)
	
	###Plot the convergence lines
	plt.axhline(0.3139, c='black')
	plt.plot(list_N, [S[0] for S in list_S_random], c='red')
	plt.plot(list_N, [S[0] for S in list_S_sobol], c='teal')
	plt.xscale('log')
	plt.show()
	plt.axhline(0.4424, c='black')
	plt.plot(list_N, [S[1] for S in list_S_random], c='red')
	plt.plot(list_N, [S[1] for S in list_S_sobol], c='teal')
	plt.xscale('log')
	plt.show()
	plt.axhline(0.0, c='black')
	plt.plot(list_N, [S[2] for S in list_S_random], c='red')
	plt.plot(list_N, [S[2] for S in list_S_sobol], c='teal')
	plt.xscale('log')
	plt.show()
	
	"""
	S, ST, model_evals = saltelli_indices("ishigami", var_names, doPrint=True)
	
	sample_sizes.append(model_evals)
	Ss.append(S)
	
	__print_ishigami_indices()
	plt.plot(sample_sizes, [S[0] for S in Ss],c='orange')
	plt.plot(sample_sizes, [S[1] for S in Ss],c='blue')
	plt.plot(sample_sizes, [S[2] for S in Ss],c='green')
	plt.axhline(0.31390519114781146, c='orange')
	plt.axhline(0.4424111447900409, c='blue')
	plt.axhline(0.0, c='green')
	plt.show()
	"""