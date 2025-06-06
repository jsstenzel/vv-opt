import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import re
import csv

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

def learn_gp_prior(wave_pts, thru_pts, err_pts=[], doPlot=False, doPrint=False, careful=False):
	X = np.array(wave_pts)
	Y = np.array(t_to_u(thru_pts)) #convert t to u-domain
	Xmin = min(wave_pts)
	Xmax = max(wave_pts)
	kernel = ConstantKernel(0.01, constant_value_bounds=(1e-4, 100.0)) * RBF(100.0, length_scale_bounds=(1.0,1000.0))
	
	n_runs = 3 if careful else 0
	
	if err_pts:
		#variances = np.array([s**2 for s in err_pts]) #I know this isnt right
		
		#convert standard deviation to u-range
		#I need to do some cleverness here. I need to convert the standard deviation from the t domain to the u domain, but that conversion has to take into account where in u-space we are.
		#So, for each std_i, make it t_to_u(t_i + std_i) - (u_i)
		stddevs_over = [t_to_u(t+s) - t_to_u(t) for t,s in zip(thru_pts,err_pts)]
		stddevs_under = [t_to_u(t) - t_to_u(t-s) for t,s in zip(thru_pts,err_pts)]
		variances = np.array([ov*un for ov,un in zip(stddevs_over,stddevs_under)])
		gpr = GaussianProcessRegressor(kernel=kernel, alpha=variances, normalize_y=True, n_restarts_optimizer=n_runs)
	else:
		gpr = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=n_runs)
		
	gpr.fit([[x] for x in X], Y) #fitting in u-range
	xplot = np.arange(Xmin,Xmax,1.0)
	mu_u, std_u = gpr.predict([[x] for x in xplot], return_std=True)
	
	#For the prior points, use the measurement locations that were provided, sorted and without diplicates
	#This makes it so that future samples from this GP have the same amount of uncertainty coming from sample locations... I think that matters
	prior_pts = sorted(list(set(wave_pts)))
	
	#Define the mean fn, in u-domain
	mean_fn_from_training = [make_mean_fn(t, xplot, mu_u) for t in prior_pts]
		
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
		plt.ylabel('throughput (u-domain)')
		plt.title('Provided data to fit')
		plt.show()
		plt.clf()
	
		#mu = u_to_t(mu_u) #need to plot in t-domain
		#std = u_to_t(std_u) #need to plot in t-domain
		plt.fill_between(xplot, (mu_u - std_u), (mu_u + std_u), facecolor='grey')
		plt.plot(xplot, mu_u, c='r')
		if err_pts:
			plt.errorbar(X, Y, yerr=err_pts, ls='none')
			
		plt.xlabel('wavelength')
		plt.ylabel('throughput (u-domain)')
		plt.title('GP fit')
		plt.show()
		plt.clf()
		"""
		
		#here's the real moneymaker
		fine_plot = np.arange(Xmin,Xmax,0.1)
		gp_std = np.array([np.sqrt(expquad_variance) for _ in fine_plot]) #plotting in u-domain
		mean_pts = np.array([make_mean_fn(f, xplot, mu_u) for f in fine_plot]) #plotting in u-domain
		plt.fill_between(fine_plot, (mean_pts - gp_std), (mean_pts + gp_std), facecolor='grey')
		plt.plot(fine_plot, mean_pts)
		plt.xlabel('wavelength')
		plt.ylabel('throughput (u-domain)')
		plt.title('The GP prior we learned')
		plt.show()
		plt.clf()
	
	return expquad_variance, expquad_ls, prior_pts, mean_fn_from_training


def learn_gp_prior_from_files(filenames, meas_std, doPlot=False, doPrint=False, careful=False):
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
	errs = [meas_std for _ in wavelengths] if meas_std > 0 else []
	variance, ls, prior_pts, mean_fn = learn_gp_prior(wavelengths, throughputs, errs, doPlot=doPlot, doPrint=doPrint)
	
	if doPlot:
		gp_prior = GaussianProcessDist1D(variance, ls, prior_pts, mean_fn)
		sample1 = gp_prior.sample()
		sample1.plot_prior(showPlot=False)
		sample2 = gp_prior.sample()
		sample2.plot_prior(showPlot=False)
		sample3 = gp_prior.sample()
		sample3.plot_prior(showPlot=False)
		plt.show()
	
	return variance, ls, prior_pts, mean_fn

def learn_save_gp_prior_to_file(save_file, filenames, meas_std, save=True, doPlot=False, doPrint=False, careful=False):
	variance, ls, prior_pts, mean_fn = learn_gp_prior_from_files(filenames, meas_std, doPlot=doPlot, doPrint=doPrint, careful=careful)
	
	#grab the mean fn points straight from the fn, so still in the u-domain
	mean_fn_pts_u = mean_fn
	
	save_file_name = save_file if save_file.endswith('.csv') else save_file+'.csv'
	if save:
		with open(save_file_name, 'w+', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow([variance])
			writer.writerow([ls])
			for p,m in zip(prior_pts, mean_fn_pts_u):
				writer.writerow([p, m])

def learn_load_gp_prior_from_file(load_file, returnObj=False, doPlot=False):
	prior_pts = []
	mean_fn_pts_u = []

	load_file_name = load_file if load_file.endswith('.csv') else load_file+'.csv'
	if not os.path.isfile(load_file_name):
		print("Couldn't find GP prior save file",load_file_name)
		sys.exit()
	
	variance = 0
	ls = 0	
	with open(load_file_name, 'r', newline='') as csvfile:
		reader = csv.reader(csvfile)
		for i,row in enumerate(reader):
			if i==0:
				variance=float(row[0])
			elif i==1:
				ls=float(row[0])
			else:
				prior_pts.append(float(row[0]))
				mean_fn_pts_u.append(float(row[1]))
	
	mean_fn_load_gp = mean_fn_pts_u
		
	if doPlot:
		fine_plot = np.arange(min(prior_pts),max(prior_pts),0.1)
		gp_std = np.array([np.sqrt(variance) for _ in fine_plot]) #plotting in u-domain
		mean_pts = np.array([make_mean_fn(f, prior_pts, mean_fn_pts_u) for f in fine_plot]) #plotting in u-domain
		plt.fill_between(fine_plot, (mean_pts - gp_std), (mean_pts + gp_std), facecolor='grey')
		plt.plot(fine_plot, mean_pts)
		plt.xlabel('wavelength')
		plt.ylabel('throughput (u-domain)')
		plt.title('GP prior '+load_file)
		plt.show()
		plt.clf()
	
	return variance, ls, prior_pts, mean_fn_load_gp


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