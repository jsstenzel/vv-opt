import sys
import os
import scipy.stats
import numpy as np

sys.path.append('../..')
#from obed.obed import *

class ProblemDefinition:
	def __init__(self, _dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors):
		self._internal_eta = _eta #eta(theta, d, x)
		self._internal_H = _H #H(theta, x)
		self._internal_G = _G #G(d, x)
		
		self.x_default = _x_default
		self.dim_d = _dim_d
		self.dim_y = _dim_y
		self.dim_theta = _dim_theta
		self.dim_x = _dim_x
		
		#priors is a list of defintions of prior distributions on the vector of thetas
		#so it's a list of pairs, first is type and second is the list of corresponding parameters
		#type checks:
		if len(_priors) != self.dim_theta:
			raise ValueError('Incorrect number of dimensions for prior distribution definition; got '+str(len(_priors))+' and expected '+str(self.dim_theta))
		for prior in _priors:
			type = prior[0]
			params = prior[1]
			
			if type not in self._allowable_prior_types.keys():
				raise ValueError("Incorrect prior probability function definition: "+str(type)+" not recognized.")
			if len(params) != self._allowable_prior_types[type]:
				raise ValueError('Wrong number of arguments for prior type '+str(type)+'; got '+str(len(params))+' and expected '+str(self._allowable_prior_types[type]))
			for param in params:
				if not isinstance(param, (int, float)):
					raise ValueError('Wrong kind of arguments for prior type '+str(type)+', need int or float.')

		#if you pass all of that,
		self.priors = _priors
	
	#optional, for documentation:
	theta_names=[]
	y_names=[]
	d_names=[]
	x_names=[]
	
	def eta(self, theta, d, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition eta: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition eta: d size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		return _internal_eta(theta, d, x)
	
	def G(self, theta, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition eta: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		return _internal_G(theta, x)
	
	def H(self, d, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition eta: d size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		return _internal_H(d, x)
	
	
	_allowable_prior_types = {'gaussian': 2, #mu, sigma
							  'uniform': 2} #left bound, right bound
	
	def prior_rvs(self, num_vals):
		vals = [] #a list length num_vals of random numbers of size dim_theta
		for prior in self.priors: ###iterate over dim_theta
			type = prior[0]
			params = prior[1]
			thetas_i = [] #wait. i might have to do this weird.
	
			#generate the rvs for this one particular theta
			if type == 'gaussian':
				mu = params[0]
				sigma = params[1]
				thetas_i = scipy.stats.norm.rvs(size=num_vals, loc=mu, scale=sigma)
			elif type == 'uniform':
				left = params[0]
				right = params[1]
				thetas_i = scipy.stats.uniform.rvs(size=num_vals, loc=left, scale=right-left) #dist is [loc, loc + scale]
			else:
				raise ValueError("prior_rvs did not expect prior type "+str(type))
				
			vals.append(thetas_i)
	
		#turn from list of rvs at each prior, to list of theta rvs
		vals = np.transpose(vals)
	
		#return
		if len(vals) == 1:
			return vals[0] #this is a list (dim_theta)
		else:
			return vals #this is a list (num_vals) of random variables (dim_theta)
	
	def prior_pdf_unnorm(self, theta):
		#evaluate and return a list of probabilities (dim_theta) for each distribution component
		probabilities = []
		
		for i,prior in enumerate(self.priors): #prior i corresponds to theta i
			type = prior[0]
			params = prior[1]
			theta_i = theta[i]
		
			if type == 'gaussian':
				mu = params[0]
				sigma = params[1]
				prob_i = scipy.stats.norm.pdf(theta_i, loc=mu, scale=sigma)
			elif type == 'uniform':
				left = params[0]
				right = params[1]
				prob_i = scipy.stats.uniform.pdf(theta_i, loc=left, scale=right-left)
			else:
				raise ValueError("prior_pdf did not expect prior type "+str(type))
				
			probabilities.append(prob_i)
			
		#need to renormalize the pdf right?
		#but unnormalized may be ok for MCMC?
		return probabilities


if __name__ == "__main__":
    #unit testing
	def eta(_theta, _d, _x):
		y = [0,0]
		y[0] =  d_[0]*_theta[0] + _theta[1] + scipy.stats.norm.rvs(loc=0, scale=_d[1], size=1)[0]
		y[1] = _x[0]*_theta[0]/_theta[1] + scipy.stats.norm.rvs(loc=0, scale=_d[0], size=1)[0]
		return y
		
	def H(_theta, _x):
		return _x[1]*_theta[0] + _theta[1]**2
		
	def Gamma(_d, _x):
		return _x[2]/(_d[0]+_d[1])
		
	toy_x_default = [1, 10, 10]

	toy_priors = [ ["uniform", [-2,-1]],
				   ["uniform", [1, 21]] ]
				   
	#maybe I can wrap prior def, default, and variable name all into one dict or list object?
	
	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	toy = ProblemDefinition(2, 2, 2, 3, eta, H, Gamma, toy_x_default, toy_priors)
	print(toy.prior_pdf_unnorm([-1.5,1.5]))
	print(toy.prior_rvs(5))
	print(toy.prior_rvs(1))
