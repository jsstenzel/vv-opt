import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import re

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF

sys.path.append('..')
from approx.gaussian_process import *



#######################################
# Machine learning for GaussianProcess
#######################################

#Note that this function is *not* about fitting a Gaussian process prior to datapoints
#Instead, it's about fitting the hyperparameters of a Gaussian process exponentional-quadratic kernel to multiple sample sets of the GP.
"""
wave_pts = all of the wavelengths of all samples, including from multiple sources
thru_pts = all of the throughputs of all samples, including from multiple sources
err_pts = all of the standard deviations of the uncertainties in Y of each sample
"""

def learn_gp_prior(wave_pts, thru_pts, err_pts=[], doPlot=False, doPrint=False):
	X = np.array(wave_pts)
	Y = np.array(thru_pts)
	Xmin = min(wave_pts)
	Xmax = max(wave_pts)
	kernel = ConstantKernel(0.0001, constant_value_bounds=(1e-10, 1.0)) * RBF(100.0, length_scale_bounds=(10.0,1000.0))
	
	if err_pts:
		variances = np.array([s**2 for s in err_pts])
		gpr = GaussianProcessRegressor(kernel=kernel, random_state=0, alpha=variances, normalize_y=False)
	else:
		gpr = GaussianProcessRegressor(kernel=kernel, random_state=0, normalize_y=False)
		
	gpr.fit([[x] for x in X], Y)
	xplot = np.arange(Xmin,Xmax,1.0)
	mu, std = gpr.predict([[x] for x in xplot], return_std=True)
	
	#For the prior points, use the measurement locations that were provided, sorted and without diplicates
	#This makes it so that future samples from this GP have the same amount of uncertainty coming from sample locations... I think that matters
	prior_pts = sorted(list(set(wave_pts)))
	
	#Define the mean fn
	def mean_fn_from_training(t):
		try:
			val = np.interp(t, xplot, mu)
		except ValueError:
			return 0.0
		return val
		
	#Interrogate the kernels for ls,variance
	params = gpr.kernel_.get_params()
	expquad_variance = params['k1'].constant_value
	expquad_ls = params['k2'].length_scale
	if doPrint:
		print(params['k1'], params['k2'], flush=True)
	
	if doPlot:	
		"""
		#Relying a little on https://stackoverflow.com/questions/74151442/how-to-incorporate-individual-measurement-uncertainties-into-gaussian-process
		plt.scatter(X, Y)
		plt.xlabel('wavelength')
		plt.ylabel('throughput')
		plt.title('Provided data to fit')
		plt.show()
		plt.clf()
	
		#plt.scatter(X, Y)
		plt.fill_between(xplot, (mu - std), (mu + std), facecolor='grey')
		plt.plot(xplot, mu, c='r')
		if err_pts:
			plt.errorbar(X, Y, yerr=err_pts, ls='none')
			
		plt.xlabel('wavelength')
		plt.ylabel('throughput')
		plt.title('GP fit')
		plt.show()
		plt.clf()
		"""
		
		#here's the real moneymaker
		fine_plot = np.arange(Xmin,Xmax,0.1)
		gp_std = np.array([np.sqrt(expquad_variance) for _ in fine_plot])
		mean_pts = np.array([mean_fn_from_training(f) for f in fine_plot])
		plt.fill_between(fine_plot, (mean_pts - gp_std), (mean_pts + gp_std), facecolor='grey')
		plt.plot(fine_plot, mean_pts)
		plt.xlabel('wavelength')
		plt.ylabel('throughput')
		plt.title('The GP prior we learned')
		plt.show()
		plt.clf()
	
	return expquad_variance, expquad_ls, prior_pts, mean_fn_from_training


def learn_gp_prior_from_files(filenames, meas_std, doPlot=False, doPrint=False):
	wavelengths = []
	throughputs = []
	
	if doPrint:
		print(filenames,flush=True)
	
	for filename in filenames:
		with open(filename, "r") as f:
			for line in f:
				words = line.split()
				if len(words)==0:
					continue
				if re.search('[a-zA-Z]', words[0]):
					continue
				if not re.search('[0-9]', words[0]):
					continue
				wavelengths.append(float(words[0].replace('#','')))
				throughputs.append(float(words[1]))
	
	#0.1% T number from Kupinski & Macleod
	errs = [meas_std**2 for _ in wavelengths] if meas_std > 0 else []
	variance, ls, prior_pts, mean_fn = learn_gp_prior(wavelengths, throughputs, errs, doPlot=doPlot, doPrint=doPrint)
	
	if doPlot:
		sample = sample_gp_prior(variance, ls, prior_pts, mean_fn)
		sample.plot_prior()
	
	return variance, ls, prior_pts, mean_fn
	

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', metavar='string', nargs='+', help='Throughput files to pull from')
	args = parser.parse_args()
	
	#0.1% T number from Kupinski & Macleod
	learn_gp_prior_from_files(args.f, 0.0001, True, True)
	
	"""
	wavelengths = []
	throughputs = []
	for filename in args.f:
		with open(filename, "r") as f:
			for line in f:
				words = line.split()
				#print(words)
				if len(words)==0:
					continue
				if re.search('[a-zA-Z]', words[0]):
					continue
				if not re.search('[0-9]', words[0]):
					continue
				wavelengths.append(float(words[0].replace('#','')))
				throughputs.append(float(words[1]))
	
	#0.1% T number from Kupinski & Macleod
	errs = [0.0001**2 for _ in wavelengths]
	variance, ls, prior_pts, mean_fn = learn_gp_prior(wavelengths, throughputs, errs, doPlot=False, doPrint=True)
	
	for _ in range(5):
		sample = sample_gp_prior(variance, ls, prior_pts, mean_fn)
		sample.plot_prior()
	"""