import os
import sys
import csv

import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from itertools import islice

sys.path.insert(0, "..")
from uq.saltelli_gsa import *

#This function evaluates a set of sample files to determine how converged it is
#Methods are based on Sarrazin et al. 2016
def total_order_convergence_tests(bootstrap_size, base_name, var_names, do_subset=0, doPrint=True):
	###Get the sample set
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
		
	#Grabbing this for later
	S_full, ST_full, __ = saltelli_indices(base_name, var_names, do_subset=do_subset, doPrint=False)
	
	if do_subset == 0:
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for row in csvreader:
				Ay.append([float(elem) for elem in row])
		#This change assumes that the data files are well-behaved; use SA_datalist_health_check first
		M = len(Ay)
		pp = len(var_names)
		By = np.zeros((M,pp+1))
		Cy = []
		
		with open(base_name+'_B.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for i,row in enumerate(csvreader):
				for e,elem in enumerate(row):
					By[i][e] = float(elem)

		for p,name in enumerate(var_names):
			Ciy = np.zeros((M,pp+1))
			with open(base_name+'_C_'+name+'.csv') as csvfile:
				csvreader = csv.reader(csvfile, delimiter=',')
				for i,row in enumerate(csvreader):
					for e,elem in enumerate(row):
						Ciy[i][e] = float(elem)
			Cy.append(Ciy)
	else:
		lim = int(do_subset/2)
		###Optionally, we can analyze less than the full set of provided samples
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for row in islice(csvreader, lim):
				Ay.append([float(elem) for elem in row])
		#This change assumes that the data files are well-behaved; use SA_datalist_health_check first
		M = len(Ay)
		pp = len(var_names)
		By = np.zeros((M,pp+1))
		Cy = []
		
		with open(base_name+'_B.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for i,row in enumerate(islice(csvreader, lim)):
				By[i] = [float(elem) for elem in row]

		for p,name in enumerate(var_names):
			Ciy = np.zeros((M,pp+1))
			with open(base_name+'_C_'+name+'.csv') as csvfile:
				csvreader = csv.reader(csvfile, delimiter=',')
				for row in islice(csvreader, lim):
					for e,elem in enumerate(row):
						Ciy[i][e] = float(elem)
			Cy.append(Ciy)
		
	if doPrint:
		print("Drawing bootstrap resamples...",flush=True)
		
	indices = [] #list of length bootstrap_size runs, each of length var_names ST's
	S_indices = [] #(just doing this to get conf intervals)
	###Generate bootstrap samples of the sample set
	for ii in range(bootstrap_size):
		if doPrint:
			print("Bootstrap sample",ii+1,'/',bootstrap_size,'...',flush=True, end='\r')
		###Bootstrap to select a new A,B,Ci
		#bootstrap sample of the original sample, with replacement 
		i_samples = np.random.randint(0, len(Ay), size=len(Ay))
		Ay_resample = [Ay[i] for i in i_samples]
		By_resample = [By[i] for i in i_samples]
		Cy_resample = [[Ciy[i] for i in i_samples] for Ciy in Cy]
		
		###Process the necessary matrices
		A = [Ay_row[:-1] for Ay_row in Ay_resample] #all but last element
		yA = [Ay_row[-1] for Ay_row in Ay_resample] #only last element
		B = [By_row[:-1] for By_row in By_resample]
		yB = [By_row[-1] for By_row in By_resample]
		
		C = []
		yC = []
		for p,Ciy in enumerate(Cy_resample):
			Ci = [Ciy_row[:-1] for Ciy_row in Ciy]
			yCi = [Ciy_row[-1] for Ciy_row in Ciy]
			C.append(Ci)
			yC.append(yCi)
		
		###Evaluate saltelli_indices on each bootstrap sample for total order indices
		f02 = np.dot(np.mean(yA),np.mean(yA))  #inner product of p-length and p-length vectors
		ST = []
		S = []
		
		yAyA = np.dot(yA, yA)      #inner product of n-length and n-length vectors
		for p,name in enumerate(var_names):
			yByCi = np.dot(yB, yC[p])  #inner product of n-length and n-length vectors
			numerator = (yByCi/M - f02)
			denominator = (yAyA/M - f02)
			if denominator == 0:
				#I think this would only ever happen if the sampler chooses all the same sample for the resample. Really unlikely for realistic calculations; just skip
				continue
			STi = 1 - numerator/denominator
			ST.append(STi)
			
			yAyCi = np.dot(yA, yC[p])
			Si = (yAyCi/M - f02)/denominator
			S.append(Si)
		
		###Add those ST (length var_names) to the boostrap sample list
		indices.append(ST)
		S_indices.append(S)
	
	#indices #list of length bootstrap_size runs, each of length var_names ST's
	data_per_index = np.array(indices).T #list of length var_names, each the list of bootstrap samples
	S_data_per_index = np.array(S_indices).T
	##############################################################################
	if doPrint:
		print("Calculating convergence metrics...",flush=True)
		
	#Magic numbers
	#small_width_threshold = 0.05 #This WAG from Sarrazin et al. 2016, defines a sufficiently-small 95% confidence
	small_width_threshold = 0.01 #but i like this one better for the Ishigami test.
	screening_threshold = 0.05 #This careful WAG from Sarazzin et al. 2016, defines a sensitivity below which ST is "negligible"
	
	###Convergence of the sensitivity indices
	#Handy method for getting conf intervals of mean of data
	def conf_interval(data, conf_level):
		#returns a tuple of the low and high bounds
		return scipy.stats.t.interval(confidence=conf_level, df=len(data)-1, loc=np.mean(data), scale=scipy.stats.sem(data))
	
	#Metric is the largest 95% confidence interval out of all of the factors
	conf_intervals = [conf_interval(index_data,0.95) for i,index_data in enumerate(data_per_index)]
	widths = [interval[1] - interval[0] for interval in conf_intervals]
	indices_metric = max(widths)
	indicesConverged = indices_metric < small_width_threshold
	
	#(do it for S too, so i can get those confidence intervals)
	S_conf_intervals = [conf_interval(index_data,0.95) for i,index_data in enumerate(S_data_per_index)]
	S_widths = [interval[1] - interval[0] for interval in S_conf_intervals]
	
	###Convergence of input factor ranking
	#handy method for calculating R
	def rank_vector(data):
		#return a list where each value is replaced by its ranking, assigning 1 to the largest value
		argsort = np.argsort(data)
		ranks = np.arange(len(data),0,-1)
		rank_dict = dict(zip(argsort, ranks))
		rankings = [rank_dict[i] for i in range(len(data))]
		return rankings
	
	#define the scaled rank distance metric between two resamples
	def rho(S_j, S_k): #calculated for every jth and kth bootstrap sample
		index_scales = []
		rank_dists = []
		for i in range(len(var_names)):
			Ri_j = rank_vector(S_j)[i]#the ranking of the ith parameter in the jth sample
			Ri_k = rank_vector(S_k)[i]#the ranking of the ith parameter in the kth sample
			Si_j = S_j[i]
			Si_k = S_k[i]
			index_scales.append( max([Si_j, Si_k])**2 )
			rank_dists.append(abs(Ri_j - Ri_k))
		denominator = sum(index_scales) #summing over i
		rho_sjk = sum([r*s for r,s in zip(rank_dists,index_scales)]) / denominator
		return rho_sjk
	
	
	#Evaluate this metric for all unique pairs of j =/= k resamples
	rho_full = []
	ranking_metric_quick = 0
	rankingConverged = True
	upper_triangular = int(0.5 * len(indices) * (len(indices)-1))
	x_percentile_95 = math.ceil(0.05*upper_triangular)
	for k,sk in enumerate(indices):
		for j,sj in enumerate(indices):
			if j>k:
				if rho(sj, sk) >= 1.0:
					ranking_metric_quick += 1
				if ranking_metric_quick >= x_percentile_95: #i can check early if the 95th percentile will be larger than 1
					rankingConverged = False
					break #if it is, fails to converge, stop calculating!
	
	#ranking_metric = np.quantile(rho_full,0.95)
	#rankingConverged = ranking_metric < 1.0
	#i.e. rankings are converged when the difference in rankings is on average less than one position
	
	###Convergence of input factor screening
	#Metric is the largest 95% confidence interval out of all of the *low-sensitivity* factors
	#This means we define any ST lower than this as being low-sensitivity
	#Identify those and only grab the indices of the low-sensitivity parameters:
	data_per_index_low = []
	screened_vars = []
	for i,name in enumerate(var_names):
		if ST_full[i]<screening_threshold:
			data_per_index_low.append(data_per_index[i])
			screened_vars.append(name)
	if len(screened_vars) > 0:
		data_per_index_low = [data_per_index[i] for i,_ in enumerate(var_names) if ST_full[i]<screening_threshold]
		#with these, calculate index convergences as above:
		conf_intervals_low = [conf_interval(index_data,0.95) for i,index_data in enumerate(data_per_index_low)]
		widths_low = [interval[1] - interval[0] for interval in conf_intervals_low]
		screening_metric = max(widths_low)
		#And now we develop a summary statistic for our confidence in the convergence of that set
		screeningConverged = screening_metric < small_width_threshold
	else:
		screeningConverged = True
		
	##############################################################################	
	#print("Calculating confidence intervals for each parameter...",flush=True)
	print("*****************************************************************")
	print(base_name, "sample set has",len(Ay)*2,"samples across A and B")
	print("Var name          ",'\t',"S",'\t',"S_conf",'\t',"ST",'\t'"ST_conf",'\t')
	for i,name in enumerate(var_names):
		print(f"{name:<18}",'\t',f"{S_full[i]:.4f}",'\t',f"{S_widths[i]:.4f}"'\t',f"{ST_full[i]:.4f}",'\t',f"{widths[i]:.4f}")
	#print("*****************************************************************", flush=True)
	
	###Report and return
	#print("*****************************************************************")
	#print(base_name, "sample set has",len(Ay)*2,"samples across A and B")
	#print("Variables:",var_names)
	#print("Full-set total-order indices:",[f"{s:.4f}" for s in ST])
	print(str(int(bootstrap_size)),"bootstrap samples made")
	print("Total order index convergence metric:",f"{indices_metric:.3f}","<",small_width_threshold)
	print("Indices converged:","TRUE" if indicesConverged else "FALSE")
	print("Total order index ranking convergence metric (quick): at least",ranking_metric_quick,"/",upper_triangular," pairs are > 1")
	print('\t',"( 95th percentile is located at the top",x_percentile_95,')')
	print("Rankings converged:","TRUE" if rankingConverged else "FALSE")
	if len(screened_vars) > 0:
		print("Screened variables (below",str(screening_threshold)+'):',screened_vars)
		print("Total order index screening convergence metric:",f"{screening_metric:.3f}")
		print("Screening converged:","TRUE" if screeningConverged else "FALSE","<",small_width_threshold)
	else:
		print("No variables screened (below",str(screening_threshold)+').')
	print("*****************************************************************", flush=True)
	return indicesConverged, False, screeningConverged
	
if __name__ == '__main__':  
	###########################
	#Full demonstration
	###########################
	
	###Set up the problem
	var_names = ['x1','x2','x3']
	dists = ['uniform', 'uniform', 'uniform']
	bounds = [[-np.pi,np.pi],[-np.pi,np.pi],[-np.pi,np.pi]]
	
	total_order_convergence_tests(1000, "ishigami", var_names, do_subset=10000, doPrint=True)
	sys.exit()