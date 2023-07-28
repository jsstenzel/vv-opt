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
from SALib.analyze import sobol, fast
import scipy.optimize as optimization

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


def sobol_saltelli(var_list, y_list, problem, conf = 0.95, doSijCalc=False, doPlot=False, doPrint=False, doReport=''):
	if len(var_list) != len(y_list):
		print("ERROR: sobol_saltelli input problem: mismatching number of model evals")
		print("Length var_list:",len(var_list))
		print("Length y_list:",len(y_list))
		exit()
		
	if len(var_list[0]) != len(var_names):
		print("ERROR: sobol_saltelli input problem: mismatching number of variables")
		print("Length var_list:",len(var_list))
		print("Length var_names:",len(var_names))
		
	X = np.array(var_list)
	Y = np.array(y_list)

	dim = len(var_names)
	problem = { #this needs to inherit from problem definition
		'names': var_names,
		'num_vars': dim,
		'bounds': [[0, 1]]*dim, #what? fix
		'dists': ['norm']*dim #what? fix
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
	plt.xlabel('Parameters', fontweight ='bold', fontsize = 15)
	plt.ylabel('S_T', fontweight ='bold', fontsize = 15)
	#plt.xticks([r + barWidth for r in range(len(varnames))], varnames)
	plt.tight_layout()
	if logplot:
		plt.yscale('log')	
	plt.show()
	print(varnames)

	
	
if __name__ == '__main__':  
	var_names = ['x1','x2','x3']
	sobol_saltelli_file('testdata.csv', var_names, doSijCalc=True, doPlot=True, doPrint=True, doReport='testreport')