import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.stats

import pymc

sys.path.append('../..')

#I need a nice simple interface to deal with these general messy things
#leave uncertainty handling to the prior fns!
	
	
"""
########## Now I need to do Gaussian process
#https://www.pymc.io/projects/docs/en/stable/api/gp.html
#the prior is returned as a TensorVariable, see https://github.com/pymc-devs/pytensor/blob/e8693bdbebca0757ab11353f121eed0c9b3acf66/pytensor/tensor/variable.py#L25
"""	

def sample_gp_prior(variance, ls, prior_pts, mean_fn):
	with pymc.Model() as model:
		cov_func = variance**2 * pymc.gp.cov.ExpQuad(input_dim=1, ls=ls)
		theta_gp = pymc.gp.Latent(cov_func=cov_func)
		#this defines the prior for the theta_gp:
		prior = theta_gp.prior("prior"+str(time.time()), np.array(prior_pts)[:, None])

	sample = GaussianProcess1D(theta_gp, prior_pts, mean_fn, prior)
	return sample

#This class is intended to work with the sample_gp_prior function
#Usually, that function is used to specify & sample a prior, returning this object
#Which contains the pymc gp which is the sample, and the mean fn for easy handling
class GaussianProcess1D:
	def __init__(self, gp, prior_pts, mean_fn, prior):
		self.gp = gp
		self.prior_pts = prior_pts
		self.mean_fn = mean_fn
		self.prior = prior
		
	#Takes the GP, measurement points and noise, and returns measured values
	#This is useful for 
	def eval_gp_cond(self, meas_pts, noise):
		with pymc.Model() as model:
			posterior = self.gp.conditional("meas"+str(time.time()), Xnew=np.array(meas_pts)[:, None])
			posterior_data = posterior.eval()
		meas_data = [y + self.mean_fn(meas_pts[i]) + scipy.stats.norm.rvs(scale=noise) for i,y in enumerate(posterior_data)]
		return meas_data
	
	#I think this is true...
	#just doing a conditional at the same prior points, with no noise applied
	def evaluate(self):
		data = self.prior.eval()
		#apply mean fn?
		for i,pt in enumerate(data):
			data[i] += self.mean_fn(self.prior_pts[i])
		return data
	
	def plot_prior(self, plotMean=True, plotDataX=[], plotDataY=[]):
		data = self.evaluate()
		plt.plot(self.prior_pts, self.evaluate(),'r')
		if plotMean:
			plt.plot(self.prior_pts, [self.mean_fn(xi) for xi in prior_pts],'g-')
		if plotDataX and plotDataY:
			plt.scatter(plotDataX, plotDataY, c='b')
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