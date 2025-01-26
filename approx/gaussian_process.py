import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.interpolate
from scipy.stats import logistic

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF

"""
Gaussian process priors and samples, specifically modeling optical throughputs
Based on the concept exhibited in this gist:
https://gist.github.com/neubig/e859ef0cc1a63d1c2ea4
"""	

def _zero_mean(x):
	return 0
	
def t_to_u(t):
	#Maps [0,1] to [-inf,inf]
	return scipy.stats.logistic.ppf(t, scale=2)
	
def u_to_t(u):
	#Maps [-inf,inf] to [0,1]
	return scipy.stats.logistic.cdf(u, scale=2)


#######################################
# Class definition
#######################################

#This class defines a Gaussian Process, such as a prior, and provides
#a method to sample from that prior to get functionals int he form of GaussianProcessSample1D
class GaussianProcessDist1D:
	def __init__(self, variance, ls, prior_pts, mean_fn=_zero_mean):
		#Variance: scaling parameter of the kernel
		self.variance = variance
		#Length scale: parameter of the kernel
		self.ls = ls
		#Mean funciton: A function defined over the whole domain, but scaled over the u-range
		self.mean_fn = mean_fn
		
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
		GP_mean = [self.mean_fn(pt) for pt in self.fine_domain]
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
		mean_fn_pts_u = [self.mean_fn(pt) for pt in self.prior_pts]
		
		save_file_name = save_file if save_file.endswith('.csv') else save_file+'.csv'
		if save:
			with open(save_file_name, 'w+', newline='') as csvfile:
				writer = csv.writer(csvfile)
				writer.writerow([self.variance])
				writer.writerow([self.ls])
				for p,m in zip(self.prior_pts, mean_fn_pts_u):
					writer.writerow([p, m])
					
	def calculate_posterior(self, meas_x, meas_y, meas_err):
		0
		#If I cared to, I could use sklearn's GaussianProcessRegression to:
		#1. Take in measurement points, values, and errors
		#2. Use fit(...) to find the new ls and var
		#3. Use predict(...) to calculate the mean function
		#4. Put these together to form a posterior GP (passing forward the prior_pts and meas_pts together as the new prior_pts)
		#5. Make a new GaussianProcessDist1D representing that posterior
		#However, I don't care to do this right now, because it isn't needed for any of llamas_snr_full's experiments

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

	def mean_fn_from_file(t):
		try:
			#scipy
			#f = scipy.interpolate.interp1d(np.array(prior_pts), np.array(mean_pts), kind=order)
			#val = f(t)
			#numpy
			val = np.interp(t, prior_pts, mean_pts_u)
		except ValueError:
			return 0.0
		return val

	if doPlot:
		plot_pts = np.array(np.linspace(min(prior_pts), max(prior_pts), len(prior_pts)*10))
		Yfit = [mean_fn_from_file(x) for x in plot_pts]
	
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
def load_gp_prior_from_file(load_file, returnObj=False, interp="linear"):
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
			
	#if interp == "linear":		#could implement other versions later
	def mean_fn_load_gp(t):
		try:
			val = np.interp(t, prior_pts, mean_fn_pts_u)
		except ValueError:
			return 0.0
		return val
	
	if returnObj:
		return GaussianProcessPrior1D(variance, ls, prior_pts, mean_fn_load_gp)
	else:
		variance, ls, prior_pts, mean_fn_load_gp

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
	eval_pts_fine = u_to_t([mean_fn(w) for w in fine_grid])

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
	mean_fn=sinabs
	
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