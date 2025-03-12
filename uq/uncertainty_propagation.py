#import argparse
import os
import sys
import shutil
import csv
import fileinput
sys.path.insert(0, "..")

import matplotlib.pyplot as plt

import numpy as np
import statistics
from scipy.stats import normaltest, norm

#datafile is a csv file where each line has format [x0, ..., xi, ..., xn, Y]
def uncertainty_prop_file(datafile, doPlot=False, doPrint=False):
	y = []
	with open(datafile) as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			y.append(float(row[-1]))

	#print(y)
	return uncertainty_prop(y, doPlot, doPrint)

#y_j: list of all model evaluation results, assumed to be single number
def uncertainty_prop(y_j, xlab="Y", c='#21aad3', saveFig='', rescaled=False, vline=[], doPlot=True, doPrint=True):
	mean = statistics.mean(y_j)
	stddev = statistics.stdev(y_j, mean) #note! This is sample stddev, not population stddev. Different n-factor in front

	score, pvalue = normaltest(y_j, nan_policy='raise')
	
	if doPrint:
		print("***Uncertainty propagation results")
		print("Sample mean:",mean)
		print("Sample standard deviation:",stddev)
		print("Gaussian fit p-value:",pvalue,"(p-value > 0.05 means its normal)")
	
	if doPlot:
		uncertainty_prop_plot(y_j, xlab, c, saveFig, rescaled, vline)
	
	return mean, stddev
	
def uncertainty_prop_plot(y_j, xlab="Y", c='#21aad3', saveFig='', rescaled=False, vline=[]):	
	"""
	###remarkably, plt.hist draws a memory error if there are outliers. it tries to list out a trillion bins.
	#Fix: prune out the worst outliers as standard practice:
	if not mean:
		mean = statistics.mean(y_j)
	if not stddev:
		stddev = statistics.stdev(y_j, mean) #note! This is sample stddev, not population stddev. Different n-factor in front
	outlier_min= mean - 5.0*stddev #5 sigma is pretty conservative
	outlier_max= mean + 5.0*stddev #5 sigma is pretty conservative
	y=[yj for yj in y_j if (yj>outlier_min and yj<outlier_max)]
	if len(y) < len(y_j):
		plt.title("(Pruned "+str(len(y_j)-len(y))+" outliers)")
	"""

	try:
		if rescaled==False:
			n, bins, patches = plt.hist(x=y_j, bins='auto', color=c, alpha=1.0, rwidth=0.85)
		else:
			n, bins, patches = plt.hist(x=y_j, bins='auto', color=c, alpha=1.0, rwidth=0.85, density=True)
	except MemoryError:
		if rescaled==False:
			n, bins, patches = plt.hist(x=y_j, bins=1000, color=c, alpha=1.0, rwidth=0.85)
		else:
			n, bins, patches = plt.hist(x=y_j, bins=1000, color=c, alpha=1.0, rwidth=0.85, density=True)	
		
	plt.grid(axis='y', alpha=0.75)
	plt.ylim(ymax=np.ceil(n.max() / 10) * 10 + n.max()*0.05)
	#plt.xticks(rotation=90)
	plt.xlabel(xlab)
	plt.ylabel("Frequency (N=" + str(len(y_j)) + ")")
	
	if vline:
		for v in vline:
			plt.axvline(vline, c='k')
	
	if saveFig=='':
		plt.show()
	else:
		if saveFig[-4:] != '.png':
			saveFig += '.png'
		plt.savefig(saveFig, format='png', bbox_inches='tight')
		plt.clf()
		plt.close()
	
def uncertainty_prop_plots(Y, xlabs=None, c='#21aad3', transpose=True, saveFig=''):
	if transpose==True:
		Y = np.transpose(Y)

	if xlabs == None:
		xlabs = ["data"+str(i) for i,_ in enumerate(Y)]
		
	if saveFig == '':
		fignames = ['' for i,_ in enumerate(Y)]
	else:
		fignames = [saveFig+"_"+str(i) for i,_ in enumerate(Y)]
		
	for i,y in enumerate(Y):
		uncertainty_prop_plot(y, xlab=xlabs[i], c=c, saveFig=fignames[i])
	
if __name__ == '__main__':  
	print("uncertainty_propagation.py: Use me to analyze and plot the results of model samples!")
	uncertainty_prop_file('testdata.csv', True, True)