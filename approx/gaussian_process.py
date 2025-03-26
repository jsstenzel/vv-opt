import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.interpolate
from scipy.stats import logistic
import csv

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF

"""
Gaussian process priors and samples, specifically modeling optical throughputs
Based on the concept exhibited in this gist:
https://gist.github.com/neubig/e859ef0cc1a63d1c2ea4
"""	

	
def t_to_u(t):
	#Maps [0,1] to [-inf,inf]
	return scipy.stats.logistic.ppf(t, scale=2)
	
def u_to_t(u):
	#Maps [-inf,inf] to [0,1]
	return scipy.stats.logistic.cdf(u, scale=2)
	
#Define the mean fn, in u-domain
def make_mean_fn(t,xplot,mu_u):
	try:
		val = np.interp(t, xplot, mu_u)
	except ValueError:
		return 0.0
	return val


#######################################
# Class definition
#######################################

#This class defines a Gaussian Process, such as a prior, and provides
#a method to sample from that prior to get functionals int he form of GaussianProcessSample1D
class GaussianProcessDist1D:
	#TODO i think this method is flawed... it would be better to dispense with the mean function,
	#just pass in the gpr that actually gets trained inlearn_gp,
	#and re-fit with additional points, as shown in test_gaussian_process_updates. (would require additional carefulness with alphas, normalization)
	#... but this works... I think. Not sure if important information gets lost when throwing out the trained gpr.
	def __init__(self, variance, ls, prior_pts, mean_fn):
		#Variance: scaling parameter of the kernel
		self.variance = variance
		#Length scale: parameter of the kernel
		self.ls = ls
		#Mean funciton: Thru values defined over prior pts, but scaled over the u-range. used for mean_fn_call
		self.mean_fn = mean_fn
		
		#I'm saving this just to reuse for if we want to get a posterior, which incolves re-training
		self.prior_pts = prior_pts
		
		#A suitably fine grid over the relevant domain -- fine enough that you don't care about smaller details, but not so fine that it takes forever to sample
		#Solution: otherwise leave it up to the user, but don't let it be finer than 1 wavelength resolution
		domain_min = min(prior_pts)
		domain_max = max(prior_pts)		
		finest_domain = np.arange(domain_min, domain_max, 1.0).tolist()
		if len(finest_domain) > len(prior_pts):
			self.fine_domain = prior_pts
		else:
			self.fine_domain = finest_domain
		
		#Make the scipy GP that you can sample from
		kernel = ConstantKernel(constant_value=variance) * RBF(length_scale=ls)
		self._gpr = GaussianProcessRegressor(kernel=kernel)
		
	def mean_fn_call(self, t):
		try:
			val = np.interp(t, self.prior_pts, self.mean_fn)
		except ValueError:
			return 0.0
		return val
	
	#Defines the radial basis function / squared exponential function that serves as the kernel
	def _rbf_kernel(self, x1, x2, ls, var):
		return var * math.exp(-1 * ((x1-x2) ** 2) / (2*ls**2))
		
	#Covariance matrix of every point in x with every other point
	def _gram_matrix(self, xs):
		mat = [[self._rbf_kernel(x1,x2,self.ls,self.variance) for x2 in xs] for x1 in xs]
		return np.asarray(mat)

	#Sample the prior in order to get a GaussianProcessSample1D
	def sample(self):
		#Sample from the Gaussian Process:
		GP_mean = [self.mean_fn_call(pt) for pt in self.fine_domain]
		#GP_cov = self._gram_matrix(self.fine_domain)
		#sample_u = np.random.multivariate_normal(GP_mean, GP_cov)
		sample_u = self._gpr.sample_y(np.array(self.fine_domain).reshape(-1,1), n_samples=1).flatten()
		sample_u_mean = [u+mean for u,mean in zip(GP_mean,sample_u)]
		
		#Now, convert it back into the t-range:
		sample_t = u_to_t(sample_u_mean)
		
		#Return the functional:
		return GaussianProcessSample1D(self.fine_domain, sample_t, u_to_t(GP_mean))
		
	def save_to_file(self, save_file, filenames, meas_std):
		#grab the mean fn points straight from the fn, so still in the u-domain
		mean_fn_pts_u = [self.mean_fn_call(pt) for pt in self.prior_pts]
		
		save_file_name = save_file if save_file.endswith('.csv') else save_file+'.csv'
		if save:
			with open(save_file_name, 'w+', newline='') as csvfile:
				writer = csv.writer(csvfile)
				writer.writerow([self.variance])
				writer.writerow([self.ls])
				for p,m in zip(self.prior_pts, mean_fn_pts_u):
					writer.writerow([p, m])
					
	def calculate_posterior(self, meas_x, meas_y, meas_err, careful=True):
		#see https://scikit-learn.org/stable/auto_examples/gaussian_process/plot_gpr_prior_posterior.html
		#code is also similar to learn_gp_prior(...)
		#1. Take in measurement points, values, and errors		
		X = np.array(meas_x)
		#convert t to u-domain
		#also need to subtract the mean, because we're fitting these points with the expectation that the prior is adjusted by the mean; meas_y needs to match the y-domain
		Y = np.array([t_to_u(y)-self.mean_fn_call(x) for x,y in zip(X,Y)]) 
		
		#1.5 Make a new GPR
		kernel = _gpr.kernel_
		n_runs = 3 if careful else 0
		
		#convert standard deviation to u-range
		#I need to do some cleverness here. I need to convert the standard deviation from the t domain to the u domain, but that conversion has to take into account where in u-space we are.
		#So, for each std_i, make it t_to_u(t_i + std_i) - (u_i)
		stddevs_over = [t_to_u(t+s) - t_to_u(t) for t,s in zip(meas_y,meas_err)]
		stddevs_under = [t_to_u(t) - t_to_u(t-s) for t,s in zip(meas_y,meas_err)]
		variances = np.array([ov*un for ov,un in zip(stddevs_over,stddevs_under)])
		new_gpr = GaussianProcessRegressor(kernel=kernel, alpha=variances, normalize_y=True, n_restarts_optimizer=n_runs)
		
		#2. Use fit(...) to find the new ls and var
		#https://stackoverflow.com/questions/74151442/how-to-incorporate-individual-measurement-uncertainties-into-gaussian-process
		new_gpr.fit([[x] for x in X], Y) #fitting in u-range
		
		#3. Use predict(...) to calculate the mean function
		#this is unfortunately necessary, because I make a new gpr object whenever I make a new GaussianProcessDist1D
		mu_u, std_u = gpr.predict([[x] for x in self.fine_domain], return_std=True)
		
		#4. Put these together to form a posterior GP
		params = gpr.kernel_.get_params()
		expquad_variance = params['k1'].constant_value
		expquad_ls = params['k2'].length_scale
		#passing forward the prior_pts and meas_pts together as the new prior_pts
		posterior_pts = sorted(list(set(list(self.prior_pts) + list(meas_x))))
		
		#Define the mean fn, in u-domain
		mean_fn_from_inference = [make_mean_fn(t, self.fine_domain, mu_u) for t in posterior_pts]
		
		#5. Make a new GaussianProcessDist1D representing that posterior
		return GaussianProcessDist1D(expquad_variance, expquad_ls, posterior_pts, mean_fn=mean_fn_from_inference)


#This class describes a functional sampled from a GP prior described by GaussianProcessDist1D
#This is not something we can calculate a posterior from; it's just a prior
class GaussianProcessSample1D:
	def __init__(self, fine_domain, thru_pts, mean_pts=[]):
		#A suitably fine grid over the relevant domain -- fine enough that you don't care about smaller details
		#Note that any information between these grid points has been lost by this point. That's why it's important to choose a good grid in the prior.
		self.fine_domain = fine_domain
		#The throughput values of the functional, defined on fine_domain
		self.thru_pts = thru_pts
		#The mean values of the GP prior, defined on fine_domain. Just for plotting/reference
		self.mean_pts = mean_pts

	#The sample of the GP prior
	#With no noise, this just evaluates the GP sample at the provided prior points
	def evaluate(self):
		return self.thru_pts
		
	#Linear interpolation is ok, because we're already (assumedly) using a fine grid for the domain
	def _interp(self, t):
		try:
			val = np.interp(t, self.fine_domain, self.thru_pts)
		except ValueError:
			return 0.0
		return val
		
	#Measurement of the sample
	#With measurement points and noise, returns measured values
	def measure(self, meas_pts, noise):
		raw_data = self._interp(meas_pts)
		if noise > 0:
			meas_data = [pt + scipy.stats.norm.rvs(scale=noise) for pt in raw_data]
		else:
			meas_data = raw_data
		return meas_data
	
	#Plot the evaluation given by evaluate()
	def plot_prior(self, plotMean=True, plotDataX=[], plotDataY=[], showPlot=True):
		plt.plot(self.fine_domain, self.thru_pts, c=np.random.rand(3))
		title="Functional"
		if plotMean and list(self.mean_pts):
			plt.plot(self.fine_domain, self.mean_pts,'k-')
			title="Sample of the Gaussian Process vs. the mean"
		if list(plotDataX) and list(plotDataY):
			plt.scatter(plotDataX, plotDataY, c='b')
		plt.title(title)
		if showPlot:
			plt.show()


##############################################################################
# Helper fns to build GaussianProcessDist1D from files
##############################################################################

#Extract prior_pts and mean_fn from a throughput file
def get_ppts_meanfn_file(filename, order=3, doPlot=False):
	prior_pts = []
	mean_pts = []
	with open(filename, "r") as f:
		for line in f:
			words = line.split()
			if len(words)==0:
				continue
			if "#" not in words[0]:
				prior_pts.append(float(words[0]))
				mean_pts.append(float(words[1]))

	#Convert to define mean function in u
	mean_pts_u = t_to_u(mean_pts)
	
	mean_fn_from_file = mean_pts_u

	if doPlot:
		plot_pts = np.array(np.linspace(min(prior_pts), max(prior_pts), len(prior_pts)*10))
		Yfit = [make_mean_fn(x, prior_pts, mean_pts_u) for x in plot_pts]
	
		plt.plot(prior_pts, mean_pts, c='k')
		plt.plot(plot_pts, Yfit, c='orange')
		plt.xlim(min(prior_pts), max(prior_pts))
		plt.title(filename)
		plt.show()

	return prior_pts, mean_fn_from_file

""" TODO right now I think these aren't or shouldn't be used, check on that
def sample_gp_from_file(filename, variance, ls, order=3):
	prior_pts, mean_fn_from_file = get_ppts_meanfn_file(filename, order)

	return sample_gp_prior(variance, ls, prior_pts, mean_fn_from_file)
	
def define_functional_from_file(filename, order=3):
	prior_pts, mean_fn_from_file = get_ppts_meanfn_file(filename, order)

	return define_functional_mean(prior_pts, mean_fn_from_file)
"""

##############################################################################
# Constructors for GaussianProcessDist1D and GaussianProcessSample1D
##############################################################################

#This is the inverse of GaussianProcessDist1D's save_to_file()
def load_gp_prior_from_file(load_file, returnObj=False, order="linear"):
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
	
	if returnObj:
		return GaussianProcessPrior1D(variance, ls, prior_pts, mean_fn_load_gp)
	else:
		return variance, ls, prior_pts, mean_fn_load_gp

#Constructor for GaussianProcessSample1D, allowing for finer polynomial interpolation
#From raw points
def define_functional(prior_pts, eval_pts, order=1):
	#Define a fine grid to interpolate onto, using half-nm wavelength steps
	fine_grid = np.arange(min(prior_pts), max(prior_pts), 0.5)
		
	#evaluate, interpolating with spline of appropriate order
	if order==1:
		eval_pts_fine = np.interp(fine_grid, prior_pts, eval_pts)
	else:
		spl = scipy.interpolate.make_interp_spline(prior_pts,eval_pts,order)
		eval_pts_fine = spl(fine_grid)

	functional = GaussianProcessSample1D(fine_grid, eval_pts_fine)
	#functional.plot_prior()
	return functional

#Constructor for GaussianProcessSample1D, allowing for finer polynomial interpolation
#From prior points (just to get the domain) and a u-range mean_fn on that domain
def define_functional_from_meanfn(prior_pts, mean_fn):
	#Define a fine grid to interpolate onto, using half-nm wavelength steps
	fine_grid = np.arange(min(prior_pts), max(prior_pts), 0.5)
		
	#evaluate, interpolating with spline of appropriate order
	eval_pts_fine = u_to_t([make_mean_fn(w,prior_pts,mean_fn) for w in fine_grid])

	functional = GaussianProcessSample1D(fine_grid, eval_pts_fine)
	#functional.plot_prior()
	return functional

#Constructor for GaussianProcessSample1D, allowing for finer polynomial interpolation
#From raw file
def define_functional_from_file(filename, order=1):
	prior_pts = []
	mean_pts = []
	with open(filename, "r") as f:
		for line in f:
			words = line.split()
			if len(words)==0:
				continue
			if "#" not in words[0]:
				prior_pts.append(float(words[0]))
				mean_pts.append(float(words[1]))
	
	return define_functional(prior_pts, mean_pts, order)

##############################################################################
# Unit test
##############################################################################

if __name__ == "__main__":	
	#define prior. Mean_fn must be on the range 0..1, and then transformed to u domain
	def sinabs(x):
		t_fn = 0.1 + 0.8*abs(math.sin(0.5*math.pi*x))
		return t_to_u(t_fn)
	
	variance=1.0 #defined on u range
	ls=.0001
	prior_pts=np.linspace(0,7,1000)
	mean_fn=[sinabs(xi) for xi in prior_pts]
	
	#print([sinabs(xi) for xi in prior_pts])

	#sample from the prior
	gp_prior = GaussianProcessDist1D(variance, ls, prior_pts, mean_fn)
	sample = gp_prior.sample()
	#print(sample.evaluate())
	
	#measure the prior
	meas_points = [0,1,2,3,4,5,6]
	noise = 0.01
	ymeas = sample.measure(meas_points, noise)
	print(ymeas)
	sample.plot_prior(plotDataX=meas_points, plotDataY=ymeas)