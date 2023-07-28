#import argparse
import os
import sys
import shutil
import csv
import fileinput

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import multiprocessing as mp
import scipy.optimize as optimization
import scipy.stats

from SALib.sample import saltelli
from SALib.analyze import sobol, fast

sys.path.insert(0, "..")
from problems.problem_definition import *

num_workers = mp.cpu_count() 


def sobol_saltelli_file(datafile, var_names, doSijCalc=False, doPlot=False, doPrint=False, doReport=''):
	y = []
	x = []
	with open(datafile) as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			x.append([float(xi) for xi in row[0:-1]])
			y.append(float(row[-1]))
	#print(x)
	#print(y)

	sobol_saltelli(x, var_names, y, doSijCalc, doPlot, doPrint, doReport)

#based on https://salib.readthedocs.io/en/latest/user_guide/basics_with_interface.html
#Generates samples, generates Y values, runs the SA
def sobol_saltelli(function, N, var_names, var_dists, var_bounds, conf = 0.95, doSijCalc=False, doPlot=False, doPrint=False, doReport=''):
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
		elif dist == 'unif' or dist == 'uniform':
			dists[j] = 'unif'
		elif dist == 'triang' or dist == 'triangle':
			dists[j] = 'triang'
		elif dist == 'norm' or dist == 'normal':
			dists[j] = 'norm'
		elif dist == 'lognorm' or dist == 'lognormal':
			dists[j] = 'lognorm'
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
	
	###generate param values
	param_values = saltelli.sample(problem, N, calc_second_order=doSijCalc)
	
	#apply the distributions unsupported by SALib
	#under advice by https://stackoverflow.com/questions/45280278/could-salib-support-other-probability-distribution-when-inputing-parameters-in-s
	for j,dist in enumerate(var_dists):
		if dist == 'gamma_mv':
			mean = var_bounds[j][0]
			variance = bounds[j][1]
			alpha = mean**2 / variance
			beta = mean / variance
			param_values[:,j] = scipy.stats.gamma.ppf(param_values[:,j], a=alpha, scale=1.0/beta)
		if dist == 'gamma_ab':
			alpha = var_bounds[j][0]
			beta = bounds[j][1]
			param_values[:,j] = scipy.stats.gamma.ppf(param_values[:,j], a=alpha, scale=1.0/beta)
	
	###run model
	Y = np.array([function(x) for x in param_values])
	
	###calculate indices
	Si = sobol.analyze(problem, Y, calc_second_order=doSijCalc, conf_level=conf, print_to_console=doPrint, parallel=False, n_processors=None) 
	if doPrint:
		print("Confidence levels on each parameter calculated at ",conf,flush=True)
		#print("S1:\t",Si['S1'],flush=True)
		#print("S1_conf:\t",Si['S1_conf'],flush=True)
		#print("ST:\t",Si['ST'],flush=True)
		#print("ST_conf:\t",Si['ST_conf'],flush=True)
	
	if doPlot:
		plot_gsa(var_names, Si=Si, logplot=True)
	
	if doReport != '':
		with open(doReport+'.csv', 'w+', newline='') as csvfile:
			writer = csv.writer(csvfile)
			for i,_ in enumerate(Si['S1']):
				row = [Si['S1'][i], Si['S1_conf'][i], Si['ST'][i], Si['ST_conf'][i]]
				if doSijCalc:
					row.append([Si['S2'][i], Si['S2_conf'][i]])
				writer.writerow(row)
				
	return Si


#This function crucially assumes that the Ys you put in were sampled according to the sobol method
def sobol_saltelli_presampled(y_list, var_names, conf = 0.95, doSijCalc=False, doPlot=False, doPrint=False, doReport=''):
	Y = np.array(y_list)

	dim = len(var_names)
	problem = { #this needs to inherit from problem definition?
		'names': var_names,
		'num_vars': dim,
		'bounds': [[0, 1]]*dim, #fix
		'dists': ['norm']*dim #fix
	}

	###calculate indices
	Si = sobol.analyze(problem, Y, calc_second_order=doSijCalc, conf_level=conf, print_to_console=doPrint, parallel=False, n_processors=None) 
	if doPrint:
		print("Confidence levels on each parameter calculated at ",conf,flush=True)
		#print("S1:\t",Si['S1'],flush=True)
		#print("S1_conf:\t",Si['S1_conf'],flush=True)
		#print("ST:\t",Si['ST'],flush=True)
		#print("ST_conf:\t",Si['ST_conf'],flush=True)
	
	if doPlot:
		#decide to logplot if the indices span 2 orders of magnitude
		checklist = [S for S in Si['ST'] if S >1e-12]
		ST_min = np.min(checklist)
		ST_max = np.max(checklist)
		plot_gsa(var_names, Si=Si, logplot=(ST_max/ST_min > 100))
	
	if doReport != '':
		with open(doReport+'.csv', 'w+', newline='') as csvfile:
			writer = csv.writer(csvfile)
			for i,_ in enumerate(Si['S1']):
				row = [Si['S1'][i], Si['S1_conf'][i], Si['ST'][i], Si['ST_conf'][i]]
				if doSijCalc:
					row.append([Si['S2'][i], Si['S2_conf'][i]])
				writer.writerow(row)
				
	return Si



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
	barWidth = 0.1
	#fig = plt.subplots(figsize =(10, 5))

	# Make the plot
	plt.bar(varnames, ST, color='orange', width = barWidth, edgecolor='grey',)
	
	# Adding Xticks
	plt.xlabel('Parameters', fontweight ='bold', fontsize = 10)
	plt.ylabel('S_T', fontweight ='bold', fontsize = 10)
	#plt.xticks([r + barWidth for r in range(len(varnames))], varnames)
	plt.tight_layout()
	if logplot:
		plt.yscale('log')	
		
	#add whiskers for conf interval	
	plt.show()
	

	
	
if __name__ == '__main__':  
	var_names = ['x1','x2','x3']
	sobol_saltelli_file('testdata.csv', var_names, doSijCalc=True, doPlot=True, doPrint=True, doReport='testreport')