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
import pickle

sys.path.append('..')
from inference.goal_based_inference import *

#Draw N samples and save in a consistent 1-line format to savefile
#Do it iteratively, to support parallelization and clustering
def bn_sampling(problem, savefile, N, buffer_rate=1, doPrint=False, sample_x=False):
	filename = savefile if savefile.endswith('.csv') else savefile+'.csv'
	
	data_buffer = []
	for i in range(N):
		if doPrint:
			print("Drawing Q & y samples",i,"...",flush=True)
		#Drawing samples theta, d
		theta_sample = problem.prior_rvs(1)
		d_sample = problem.sample_d(1)
		if sample_x:
			x_sample = problem.sample_x(1)
		else:
			x_sample=[]
			
		#Model propagation for Q, y
		qoi_train = problem.H(theta_sample, x_sample)
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
			elif not do_subset or len(Q) < do_subset:
				try:
					ygrab = [float(e) for e in row[:problem.dim_y]]
					dgrab = [float(e) for e in row[problem.dim_y:-1]]
					Qgrab = float(row[-1])
				except ValueError: #recently im seeing some '' values in y? hopefully this avoids that ugliness
					continue
				y.append(ygrab) #should be length dim_y
				d.append(dgrab) #should be length dim_d = row - dim_y - 1
				Q.append(Qgrab)
			else:
				break
	
	#zip it together at the end, in case i ever need them separate for something later
	yd = [y_i + d_i for y_i,d_i in zip(y,d)]

	return Q, yd

#Call bn_load_samples and gbi_train_model
def bn_train_from_file(problem, savefile, do_subset=0, ncomp=0, pca=None, doPrint=False):
	###Load file
	qoi_train, y_d_train = bn_load_samples(problem, savefile, doPrint, do_subset)
	#y_d_train = [[yd[0], yd[1]] for yd in y_d_train] #stupid cut for speed
	
	###Train model
	gmm = gbi_train_model(qoi_train, y_d_train, verbose=2, ncomp=ncomp, pca=pca, careful=True)
	
	###Print and return
	if doPrint:
		print("Trained GMM with",len(qoi_train),"samples from",savefile)
		print(gmm,len(gmm.means_[0]),"dimensions", flush=True)
		
	return gmm
	

#Save the GMM to file for easy grabbing later
def bn_save_gmm(gmm, gmm_file):
	filename = gmm_file if gmm_file.endswith('.pkl') else gmm_file+'.pkl'
	
	with open(filename, 'wb') as file:
		pickle.dump(gmm, file)
	
	print("GMM saved to",filename,flush=True)

#Load GMM from file
def bn_load_gmm(gmm_file):
	filename = gmm_file if gmm_file.endswith('.pkl') else gmm_file+'.pkl'
	
	with open(filename, 'rb') as file:
		data_gmm = pickle.load(file)

	return data_gmm

#Read and interpret the savefile
#then, return only a randomized sub-sample of y
def bn_load_y(problem, savefile, do_subset=0, doPrint=False, doDiagnostic=False):
	filename = savefile if savefile.endswith('.csv') else savefile+'.csv'

	###Make sure the files exist
	if not os.path.isfile(filename):
		print("File",filename,"is missing")
		sys.exit()
		
	###Safely read out all of the samples into matrices
	y = []
	
	if doPrint:
		print("Reading the data files...",flush=True)
	
	with open(filename) as csvfile:
		csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
		for l,row in enumerate(csvreader):
			if len(row) != problem.dim_y + problem.dim_d + 1:
				if doDiagnostic:
					print("Warning: dropped line",l+1,"(length "+str(len(row))+' expected', str(problem.dim_y + problem.dim_d + 1)+')',"from",filename)
			elif not do_subset or len(y) < do_subset:
				try:
					ygrab = [float(e) for e in row[:problem.dim_y]]
				except ValueError: #recently im seeing some '' values in y? hopefully this avoids that ugliness
					continue
				y.append(ygrab) #should be length dim_y
			else:
				break
	
	return y

#take in a list of samples
#return a subset of that list of length N_draw, with randomized selections
def bn_random_subset(samples, N_draw, allowDuplicates=True, seed=None):
	###if you don't provide a seed, make it random
	if seed:
		np.random.seed(seed)
	
	###prepare the subset index list
	if allowDuplicates:
		i_subset = np.random.randint(0, len(samples), size=N_draw)
	else:
		i_subset = random.sample(range(0, len(samples)), N_draw)
	
	subset = [samples[i] for i in i_subset]
	return subset