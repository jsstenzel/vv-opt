#import argparse
import os
import sys
import shutil
import csv
import fileinput

import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import scipy.optimize as optimization
import scipy.stats

from SALib.sample import saltelli
from SALib.analyze import sobol, fast

sys.path.insert(0, "..")
from problems.problem_definition import *

num_workers = mp.cpu_count() 


def sobol_saltelli_file(datafile, var_names, var_dists, var_bounds, doSijCalc=False, doPlot=False, doPrint=False, writeReport=''):
	y = []
	x = []
	with open(datafile) as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			x.append([float(xi) for xi in row[0:-1]])
			y.append(float(row[-1]))
	#print(x)
	#print(y)
	
	problem = sobol_preprocess(var_names, var_dists, var_bounds)

	###calculate indices
	sobol_analyze_samples(problem, y, doSijCalc, conf, doPrint, writeReport)
	
	if doPlot:
		plot_gsa(var_names, Si=Si, logplot=True)
	
	return Si


#based on https://salib.readthedocs.io/en/latest/user_guide/basics_with_interface.html
#Generates samples, generates Y values, runs the SA
def sobol_saltelli(function, N, var_names, var_dists, var_bounds, conf = 0.95, doSijCalc=False, doPlot=False, doPrint=False, writeSamples='', writeReport=''):
	#check dimensions
	#Handle the dists:
	problem = sobol_preprocess(var_names, var_dists, var_bounds)
	
	#Handle the dists:
	###generate param values
	#apply the distributions unsupported by SALib
	###run model
	#write samples
	param_values, Y = sobol_generate_samples(function, N, problem, var_dists, var_bounds, doSijCalc=False, doPrint=False, writeSamples=writeSamples)
	
	###calculate indices
	Si = sobol_analyze_samples(problem, Y, doSijCalc, conf, doPrint, writeReport)
	
	if doPlot:
		plot_gsa(var_names, Si=Si, logplot=True)
	
	return Si

#Take in the provided names/dists/bounds,
#check if they're valid
#and make the subsistutions necessary for unsupported dists
def sobol_preprocess(var_names, var_dists, var_bounds):
	#check dimensions
	if len(var_names) != len(var_dists) or len(var_names) != len(var_bounds):
		print("ERROR: sobol_saltelli input problem: mismatching number of model evals")
		print("Length var_names:",len(var_names))
		print("Length var_dists:",len(var_dists))
		print("Length var_bounds:",len(var_bounds))
		exit()

	#Handle the dists:
	dists = var_dists
	bounds = var_bounds
	for j,dist in enumerate(var_dists):
		if dist == 'gamma_mv' or dist == 'gamma_ab':
			dists[j] = 'unif'
			bounds[j] = [0,1]
		elif dist == 'beta':
			dists[j] = 'unif'
			bounds[j] = [0,1]
		elif dist == 'unif' or dist == 'uniform':
			dists[j] = 'unif'
		elif dist == 'triang' or dist == 'triangle':
			dists[j] = 'triang'
		elif dist == 'norm' or dist == 'normal' or dist == 'gaussian':
			dists[j] = 'norm'
		elif dist == 'lognorm' or dist == 'lognormal':
			dists[j] = 'lognorm'
		elif dist == 'nonrandom':
			dists[j] = 'unif'
			val = bounds[j][0]
			bounds[j] = [val,val+1e-6]
		else:
			print("sobol_saltelli dist not recognized:",dist)
			exit()

	dim = len(var_names)
	problem = { #this needs to inherit from problem definition?
		'names': var_names,
		'num_vars': dim,
		'bounds': bounds,
		'dists': dists
	}
	
	return problem


def sobol_generate_samples(function, N, problem, var_dists, var_bounds, doSijCalc=False, doPrint=False, writeSamples=''):	
	###generate param values
	if doPrint:
		print("Generating Sobol sequence sample...", flush=True)
	param_values = saltelli.sample(problem, N, calc_second_order=doSijCalc)
	print(N, len(param_values))
	
	#apply the distributions unsupported by SALib
	#under advice by https://stackoverflow.com/questions/45280278/could-salib-support-other-probability-distribution-when-inputing-parameters-in-s
	for j,dist in enumerate(var_dists):
		if dist == 'gamma_mv':
			mean = var_bounds[j][0]
			variance = var_bounds[j][1]
			alpha = mean**2 / variance
			beta = mean / variance
			param_values[:,j] = scipy.stats.gamma.ppf(param_values[:,j], a=alpha, scale=1.0/beta)
		if dist == 'gamma_ab':
			alpha = var_bounds[j][0]
			beta = var_bounds[j][1]
			param_values[:,j] = scipy.stats.gamma.ppf(param_values[:,j], a=alpha, scale=1.0/beta)
		if dist == 'beta':
			a = var_bounds[j][0]
			b = var_bounds[j][1]
			param_values[:,j] = scipy.stats.beta.ppf(param_values[:,j], a=a, b=b)
		if dist == 'nonrandom':
			param_values[:,j] = np.repeat(var_bounds[j][0], len(param_values[:,j]))

	###run model
	if doPrint:
		print("Running the model over the Sobol samples...", flush=True)
		Yget = []
		for i,x in enumerate(param_values):
			Yget.append(function(x))
			print(str(i+1)+"/"+str(len(param_values)),'\t', flush=True, end='\r')
		Y = np.array(Yget)
	else:
		Y = np.array([function(x) for x in param_values])

	#write samples
	if writeSamples != '':
		#make a csv file where each line has format [x0, ..., xi, ..., xn, Y]
		writelist = [np.append(val_list, yi) for val_list, yi in zip(param_values,Y)]
		
		with open(writeSamples, 'w+', newline='') as csvfile:
			writer = csv.writer(csvfile)
			for row in writelist:
				writer.writerow(row)
				
	return param_values, Y


def sobol_analyze_samples(problem, Y, doSijCalc, conf, doPrint=False, writeReport=''):
	###calculate indices
	if doPrint:
		print("Analyzing Sobol indices...", flush=True)
	Si = sobol.analyze(problem, Y, calc_second_order=doSijCalc, conf_level=conf, print_to_console=doPrint, parallel=False, n_processors=None) 
	if doPrint:
		print("Confidence levels on each parameter calculated at ",conf,flush=True)
		#print("S1:\t",Si['S1'],flush=True)
		#print("S1_conf:\t",Si['S1_conf'],flush=True)
		#print("ST:\t",Si['ST'],flush=True)
		#print("ST_conf:\t",Si['ST_conf'],flush=True)
	
	if writeReport != '':
		with open(writeReport, 'w+', newline='') as csvfile:
			writer = csv.writer(csvfile)
			for i,_ in enumerate(Si['S1']):
				row = [Si['S1'][i], Si['S1_conf'][i], Si['ST'][i], Si['ST_conf'][i]]
				if doSijCalc:
					row.append([Si['S2'][i], Si['S2_conf'][i]])
				writer.writerow(row)
	
	return Si



def sobol_resample(sample_datafile, function, N, problem, var_dists, var_bounds, doSijCalc=False, doPrint=False, doUpdateFile=True):	
	###Get the currently-existing samples
	y_read = []
	x_read = []
	with open(sample_datafile, 'r') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			x_read.append([float(xi) for xi in row[0:-1]])
			y_read.append(float(row[-1]))
	#print(x_read)
	#print(y_read)
	
	###Compare N_prev to N
	length_prev = len(y_read)
	
	#From documentation:
	#If calc_second_order is False, the resulting matrix has N * (D + 2) rows, where D is the number of parameters.
	#If calc_second_order is True, the resulting matrix has N * (2D + 2) rows
	d = len(x_read[0])
	N_prev = int(length_prev/(2*d+2)) if doSijCalc else int(length_prev/(d+2))
	#expected_length = N_prev*(2*d+2) if doSijCalc else N_prev*(d+2)
	
	if not (N_prev > 0 and (N_prev & (N_prev - 1)) == 0):
		print("Provided sample file",sample_datafile,"is not a power of 2: N_prev =",N_prev)
		exit()
	if N <= N_prev:
		print("N ("+str(N)+") is expected to be a larger power of 2 than N_prev ("+str(N_prev)+").")
		exit()
	
	###generate N param values, even though that involves some redundancy
	#best i can do if I stick with this library
	if doPrint:
		print("Re-generating Sobol sequence sample...", flush=True)
	param_values = saltelli.sample(problem, N, calc_second_order=doSijCalc)
	
	###apply the distributions unsupported by SALib
	#under advice by https://stackoverflow.com/questions/45280278/could-salib-support-other-probability-distribution-when-inputing-parameters-in-s
	for j,dist in enumerate(var_dists):
		if dist == 'gamma_mv':
			mean = var_bounds[j][0]
			variance = var_bounds[j][1]
			alpha = mean**2 / variance
			beta = mean / variance
			param_values[:,j] = scipy.stats.gamma.ppf(param_values[:,j], a=alpha, scale=1.0/beta)
		if dist == 'gamma_ab':
			alpha = var_bounds[j][0]
			beta = var_bounds[j][1]
			param_values[:,j] = scipy.stats.gamma.ppf(param_values[:,j], a=alpha, scale=1.0/beta)
		if dist == 'beta':
			a = var_bounds[j][0]
			b = var_bounds[j][1]
			param_values[:,j] = scipy.stats.beta.ppf(param_values[:,j], a=a, b=b)
		if dist == 'nonrandom':
			param_values[:,j] = np.repeat(var_bounds[j][0], len(param_values[:,j]))

	###run model only for the new points
	if doPrint:
		print("Running the model over the new Sobol samples...", flush=True)
	Yget = []	
	for ii,xx in enumerate(param_values):
		if ii < length_prev:
			Yget.append(y_read[ii])
			#Grab the already-generated model eval -- this is how we can save time when extending a sequence
		else:
			Yget.append(function(xx))
			if doPrint:
				print(str(ii+1)+"/"+str(len(param_values)),'\t', flush=True, end='\r')
	Y = np.array(Yget)

	###write these over the provided samples
	if doUpdateFile:
		#make a csv file where each line has format [x0, ..., xi, ..., xn, Y]
		writelist = [np.append(val_list, yi) for val_list, yi in zip(param_values,Y)]
		
		with open(sample_datafile, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			for row in writelist:
				writer.writerow(row)
				
	return param_values, Y



def plot_gsa(varnames, Si=None, filename=None, logplot=False):
	S1, S1_conf, ST, ST_conf = [],[],[],[]
	if Si != None:
		S1	  = Si['S1']
		S1_conf = Si['S1_conf']
		ST	  = Si['ST']
		ST_conf = Si['ST_conf']
	elif filename != None:
		with open(filename) as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for row in csvreader:
				S1.append(float(row[0]))
				S1_conf.append(float(row[1]))
				ST.append(float(row[2]))
				ST_conf.append(float(row[3]))
				
	
	# set width of bar
	barWidth = 0.2
	#fig = plt.subplots(figsize =(10, 5))

	# Make the plot
	plt.bar(varnames, ST, color='orange', width = barWidth, edgecolor='grey',)
	
	# Adding Xticks
	plt.xlabel('Parameters', fontweight ='bold', fontsize = 10)
	plt.ylabel('S_T', fontweight ='bold', fontsize = 10)
	#plt.xticks([r + barWidth for r in range(len(varnames))], varnames)
	plt.xticks(rotation=90)
	plt.tight_layout()
	if logplot:
		plt.yscale('log')	
		
	#add whiskers for conf interval	
	plt.errorbar(varnames, ST, yerr=ST_conf, fmt='|', color='k')
	plt.grid(axis = 'y')
	
	plt.show()

	

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

if __name__ == '__main__':  
	###########################
	#Full demonstration
	###########################
	
	#Set up the problem
	var_names = ['x1','x2','x3']
	dists = ['norm', 'norm', 'norm']
	bounds = [[-np.pi,np.pi],[-np.pi,np.pi],[-np.pi,np.pi]]
	
	dim = len(var_names)
	problem = { #this needs to inherit from problem definition?
		'names': var_names,
		'num_vars': dim,
		'bounds': bounds,
		'dists': dists
	}
	
	conf = 0.95
	sample_file = "test_samples.csv"
	plist = np.arange(1,16)
	Si_list = []
	
	###Do loop of calculating indices and re-sampling
	for j,p in enumerate(plist):
		N = 2**p
		
		#Generate/re-generate samples
		if j==0:
			###Generate the initial samples	
			param_values, Y = sobol_generate_samples(__ishigami, 2**1, problem, dists, bounds, doSijCalc=False, doPrint=True, writeSamples=sample_file)
		else:
			param_values, Y = sobol_resample(sample_file, __ishigami, N, problem, dists, bounds, doSijCalc=False, doPrint=True, doUpdateFile=True)
		
		#calculate indices
		Si = sobol_analyze_samples(problem, Y, doSijCalc=False, conf=conf, doPrint=True)
		Si_list.append(Si["S1"])
		
	
	###Plot convergence
	for i,var_name in enumerate(var_names):
		plt.plot(plist, [Si[i] for Si in Si_list])
	plt.legend(var_names)
	plt.show()