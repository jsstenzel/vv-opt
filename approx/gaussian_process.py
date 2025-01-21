import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.stats
import scipy.interpolate
from scipy.stats import logistic

import pymc

sys.path.append('../..')

#I need a nice simple interface to deal with these general messy things
#leave uncertainty handling to the prior fns!
	
	
"""
########## Now I need to do Gaussian process
#https://www.pymc.io/projects/docs/en/stable/api/gp.html
#the prior is returned as a TensorVariable, see https://github.com/pymc-devs/pytensor/blob/e8693bdbebca0757ab11353f121eed0c9b3acf66/pytensor/tensor/variable.py#L25
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
# Generators for GaussianProcess1D
#######################################

def sample_gp_prior(variance, ls, prior_pts, mean_fn=_zero_mean):
	with pymc.Model() as model:
		cov_func = variance * pymc.gp.cov.ExpQuad(input_dim=1, ls=ls)
		theta_gp = pymc.gp.Latent(cov_func=cov_func)
		#this defines the prior for the theta_gp:
		prior = theta_gp.prior("prior"+str(time.time()), np.array(prior_pts)[:, None])

	sample = GaussianProcess1D(theta_gp, prior_pts, mean_fn, prior)
	return sample
	
def define_functional(prior_pts, eval_pts, order=1):
	#Convert to define mean function in u
	eval_pts_u = t_to_u(eval_pts)
		
	def mean_fn_from_pts(t):
		try:
			val = np.interp(t, prior_pts, eval_pts_u)
		except ValueError:
			return 0.0
		return val

	sample = GaussianProcess1D(None, prior_pts, mean_fn_from_pts, None)
	return sample
	
def define_functional_mean(prior_pts, mean_fn=_zero_mean):
	sample = GaussianProcess1D(None, prior_pts, mean_fn, None)
	return sample
	
#######################################
# Helper fns to do it from files
#######################################
	
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
	
def sample_gp_from_file(filename, variance, ls, order=3):
	prior_pts, mean_fn_from_file = get_ppts_meanfn_file(filename, order)

	return sample_gp_prior(variance, ls, prior_pts, mean_fn_from_file)
	
def define_functional_from_file(filename, order=3):
	prior_pts, mean_fn_from_file = get_ppts_meanfn_file(filename, order)

	return define_functional_mean(prior_pts, mean_fn_from_file)
	
#######################################
# Class definition
#######################################

#This class is intended to work with the sample_gp_prior function
#Usually, that function is used to specify & sample a prior, returning this object
#Which contains the pymc gp which is the sample, and the mean fn for easy handling
class GaussianProcess1D:
	def __init__(self, gp, prior_pts, mean_fn, prior):
		self.gp = gp
		self.prior_pts = prior_pts
		self.mean_fn = mean_fn #this is defined in u-domain
		self.prior = prior
		
		if gp==None and prior==None:
			self.static_flag = True
		else:
			self.static_flag = False
		
	#Takes the GP, measurement points and noise, and returns measured values
	#This is useful for 
	def eval_gp_cond(self, meas_pts, noise):
		if not self.static_flag:		
			with pymc.Model() as model:
				posterior = self.gp.conditional("meas"+str(time.time()), Xnew=np.array(meas_pts)[:, None])
				posterior_data = posterior.eval()
		else:
			posterior_data = [0]*len(meas_pts)
		if noise > 0:
			meas_data = [y + self.mean_fn(meas_pts[i]) + scipy.stats.norm.rvs(scale=noise) for i,y in enumerate(posterior_data)]
		else:
			meas_data = [y + self.mean_fn(meas_pts[i]) for i,y in enumerate(posterior_data)]
		#Mean fn is innately in u, nNeed to convert back to t
		meas_data_t = u_to_t(meas_data)
		return meas_data_t
	
	#I think this is true...
	#just doing a conditional at the same prior points, with no noise applied
	def evaluate(self):
		if not self.static_flag:
			data = self.prior.eval()
		else:
			data = [0]*len(self.prior_pts)
		#apply mean fn?
		for i,pt in enumerate(data):
			data[i] += self.mean_fn(self.prior_pts[i])
		#Mean fn is innately in u, need to convert back to t
		data_t = u_to_t(data)
		return data_t
	
	def plot_prior(self, plotMean=True, plotDataX=[], plotDataY=[], showPlot=True):
		data = self.evaluate()
		plt.plot(self.prior_pts, data, c=np.random.rand(3))
		if plotMean:
			plt.plot(self.prior_pts, u_to_t([self.mean_fn(xi) for xi in self.prior_pts]),'k-')
		if list(plotDataX) and list(plotDataY):
			plt.scatter(plotDataX, plotDataY, c='b')
		plt.title("Sample of the Gaussian Process")
		if showPlot:
			plt.show()


if __name__ == "__main__":
	###Test 1
	"""
	with pymc.Model() as model:
		cov_func = pymc.gp.cov.ExpQuad(input_dim=1, ls=0.1)
		theta_gp = pymc.gp.Latent(cov_func=cov_func)
		
		X = np.linspace(0, 1, 10)[:, None]
		samples = []
		for i in range(10):
			f = theta_gp.prior("f"+str(i), X) #one sample of the GP
			samples.append(f.eval())
	"""
	
	###Test 2
	"""
	#for sample in samples:
	#	plt.plot(sample)
	#plt.show()
		
	#ok, now how do i update a prior? Do i use the Latent.conditional method?
	#Once i have this, i can do the experiment model i think...
	
	#define a prior in the shape of a parabola
	prior_y = []
	prior_x = np.linspace(0, 1, 11)[:, None]
	prior2_x = np.linspace(0, 1, 20)[:, None]
	for x in prior_x:
		prior_y.append(-(x-5)**2 + 20)
	prior_pts = [[x,p] for x,p in zip(prior_x[:],prior_y)]
	print(prior_pts)
	sigma=10
	
	#https://www.pymc.io/projects/docs/en/stable/learn/core_notebooks/Gaussian_Processes.html#additive-gp
	with pymc.Model() as model:
		cov_func = sigma**2 * pymc.gp.cov.ExpQuad(input_dim=1, ls=0.1)
		theta_gp = pymc.gp.Latent(cov_func=cov_func) #mean_func=
		prior = theta_gp.prior("prior", prior_x)
		
		x_data = np.linspace(0, 1, 101)[:, None] #np.array([.1,.15,.7,.75])[:, None]
		posterior = theta_gp.conditional("posterior", Xnew=x_data) #???
		prior2 = theta_gp.prior("prior2", prior_x)
	
	print(prior.eval())
	print(posterior.eval())
	plt.plot(prior_x, prior.eval(),'r')
	plt.plot(prior_x, prior2.eval(),'g')
	plt.scatter(x_data, posterior.eval())
	plt.show()
	#this experiment shows that calling prior newly sets the prior of the gp
	#then, when you call the conditional, its kind of like taking measurements at those points, and it gives you noisy measurements around that prior
	#so now, all i have to do is create a class that systematizes this into something that lets me define a prior and an eta,
	#and also includes the option to set any kind of mean prior as a function
		
	
	#and how do i evaluate the probability density of the GP? Is this actually necessary for doing goal-based inference?
	#nah!
	"""
	
	###Real Test
	#define prior
	def parabola(x):
		return -(x-4)**2 + 20
	
	variance=.5
	ls=0.05
	prior_pts=[0,1,2,3,4,5,6,7]
	mean_fn=parabola
	
	print([parabola(xi) for xi in prior_pts])

	#sample from the prior
	sample = sample_gp_prior(variance, ls, prior_pts, mean_fn)
	print(sample.evaluate())
	
	#measure the prior
	meas_points = [0,2,4,6]
	noise = 2.0
	ymeas = sample.eval_gp_cond(meas_points, noise)
	print(ymeas)
	sample.plot_prior(plotDataX=meas_points, plotDataY=ymeas)