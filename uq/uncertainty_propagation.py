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
	
def uncertainty_prop_plots(Y, xlabs=None, c='#21aad3', transpose=True, saveFig='', vline_per_plot=[]):
	if transpose==True:
		Y = np.transpose(Y)

	if xlabs == None:
		xlabs = ["data"+str(i) for i,_ in enumerate(Y)]
		
	if saveFig == '':
		fignames = ['' for i,_ in enumerate(Y)]
	else:
		fignames = [saveFig+"_"+str(i) for i,_ in enumerate(Y)]
		
	for i,y in enumerate(Y):
		vline = [] if not vline_per_plot else [vline_per_plot[i]]
		uncertainty_prop_plot(y, xlab=xlabs[i], c=c, saveFig=fignames[i], vline=vline)
	
#this function is intended to take in a set of n_mc results from a big Monte Carlo run
#and show you how the mean and confidence interval of that one run change over n-factor
#pulling from my past work proj1_p2_16940.py
def mc_plot_trace(samples, plotConf=True, showPlot=True, c='black', a=1.0, doLog=False, doEvery=1):
	#generate the trace statistics
	n_mc = len(samples)
	n_list = []
	means = []
	confs = []
	for i in range(n_mc):
		if i % doEvery == 0:
			print(i,flush=True,end='\r')
			ni = i+1
			if i == 0:
				0 #skip this, no sense in plotting it
			else:
				n_list.append(ni)
				means.append(np.mean(samples[:ni]))
				confs.append(np.sqrt(statistics.variance(samples[:ni]) / ni) * 1.96)
	print()
	
	#then, the convergence plot
	plt.title("Trace of MC estimator as n increases")
	plt.xlabel("n_mc")
	plt.ylabel("MC mean")
	if doLog:
		plt.xscale('log')
	plt.plot(n_list,means, c=c, alpha=a)
	#then plot the confidence interval
	if plotConf:
		plt.plot(n_list,[mean-confs[i] for i,mean in enumerate(means)], c='r')
		plt.plot(n_list,[mean+confs[i] for i,mean in enumerate(means)], c='r')
	
	if showPlot:
		plt.show()
		
	return means


#do the above, with multiple MC results, plotting them all and using them all for confidence interval
#this was helpful: 
#https://stats.stackexchange.com/questions/525063/trying-to-calculate-confidence-intervals-for-a-monte-carlo-estimate-of-pi-what
def mc_plot_trace_ensemble(main_sample, multi_samples, doLog=False, savePlot=False, doEvery=1):
	n_mc = len(multi_samples[0])
	n_runs = len(multi_samples)
	
	#first, plot the traces of all of the sample sets
	multi_means = [None]*n_runs
	for j in range(n_runs):
		print(str(j),"runs...\t", end='\r', flush=True)
		means = mc_plot_trace(multi_samples[j], plotConf=False, showPlot=False, c='gray', a=0.5, doLog=doLog, doEvery=doEvery)
		multi_means[j] = means #these means are indexed at a rate of doEvery
		#print(len(means))
		
	#get the average highest-n prediction of the MC result
	best_estimate = np.mean([means[-1] for means in multi_means])
		
	#then, plot the main sample we're given (separate from the CI calculation!)
	mc_plot_trace(main_sample, plotConf=False, showPlot=False, c='black', a=1.0, doLog=doLog, doEvery=doEvery)
		
	n_list = []
	#argh, get the n_list that matches what mc_plot_trace gives
	for i in range(n_mc):
		if i % doEvery == 0 and i != 0:
			n_list.append(i+1)
		
	#then, use all of the sample set means to find the standard error, use to calculate ci
	overall_ci = []
	for i,n in enumerate(n_list):
		sample_means_at_i = [trace[i] for trace in multi_means] #the trace value at n=1 for all experiments in the ensemble
		#print("means at",i,sample_means_at_i)
		std = np.std(sample_means_at_i, ddof=1)
		se1 = std/np.sqrt(n_runs)
		se2 = std/np.sqrt(n)
		#print("standard error",std,se1,se2)
		overall_ci.append(std * 1.96) #standard error of the means, times z
		#weirdly, n_runs is missing from that expreesion... is it already included in the variance calc here??
	
	#lastly, carefully plot those cis around the trace we plotted in black
	#TODO change this so that its around the mean of all traces
	plt.plot(n_list,[best_estimate-ci for ci in overall_ci], c='r')
	plt.plot(n_list,[best_estimate+ci for ci in overall_ci], c='r')
	if savePlot:
		plt.savefig("mc_plot_trace.png", bbox_inches='tight', transparent=True)
	else:
		plt.show()
	
	print("n:\t CI:")
	for n,ci in zip(n_list,overall_ci):
		print(n, ci)
	print("Best estimate for MC result (at n="+str(n_list[-1])+"):",best_estimate,flush=True)
	
	return best_estimate

#performs bootstrapping on one MC run, and then does mc_plot_trace_ensemble
def mc_plot_trace_bootstrap(samples, n_bootstrap, doLog=False, savePlot=False, doEvery=1):
	bootstrap_samples = [None]*n_bootstrap
	for j in range(n_bootstrap):
		###Make bootstrap sample
		i_samples = np.random.randint(0, len(samples), size=len(samples))
		bootstrap_sample = [samples[i] for i in i_samples]
		bootstrap_samples[j] = bootstrap_sample
		
	return mc_plot_trace_ensemble(samples, bootstrap_samples, doLog=doLog, savePlot=savePlot, doEvery=doEvery)
	
if __name__ == '__main__':  
	print("uncertainty_propagation.py: Use me to analyze and plot the results of model samples!")
	uncertainty_prop_file('testdata.csv', True, True)
