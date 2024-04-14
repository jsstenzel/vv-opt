import sys
import os
import scipy.stats
import scipy.interpolate
import numpy as np
import math

sys.path.append('..')
from problems.gaussian_process import *


class ProblemDefinition:
	def __init__(self, _eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs):
		self._internal_eta = _eta #eta(theta, d, x)
		self._internal_H = _H #H(theta, x)
		self._internal_G = _G #G(d, x)
		
		#get dimensions, set default parameters
		self.dim_d = len(_d_defs)
		self.dim_y = len(_y_defs)
		self.dim_theta = len(_theta_defs)
		self.dim_x = len(_x_defs)
		self.x_default = [default for name,dist,mask,default in _x_defs]
		
		#hack the x_defs to save wasted effort:
		for xdef in _x_defs:
			if xdef[1] == []:
				xdef[1] = ["nonrandom", [xdef[3]]]
		
		#priors is a list of defintions of prior distributions on the vector of thetas
		#so it's a list of pairs, first is type and second is the list of corresponding parameters
		#type checks:
		for _,prior,mask in _theta_defs + _d_defs + [x[:-1] for x in _x_defs]:
			type = prior[0]
			params = prior[1]
			
			if type not in self._allowable_prior_types.keys():
				raise ValueError("Incorrect prior probability function definition: "+str(type)+" not recognized.")
			if len(params) != self._allowable_prior_types[type]:
				raise ValueError('Wrong number of arguments for prior type '+str(type)+'; got '+str(len(params))+' and expected '+str(self._allowable_prior_types[type]))
			if not (mask in ['continuous','discrete','functional','object']):
				raise ValueError('Variable mask not an expected value: '+str(mask))

		#if you pass all of that,
		self.priors  = [prior for _,prior,_ in _theta_defs]
		self.d_dists = [dist for _,dist,_ in _d_defs]
		self.x_dists = [dist for _,dist,_,_ in _x_defs]
		
		self.theta_masks = [mask for _,_,mask in _theta_defs]
		self.d_masks     = [mask for _,_,mask in _d_defs]
		self.x_masks     = [mask for _,_,mask,_ in _x_defs]
		
		#for documentation:
		self.theta_names=[name for name,_,_ in _theta_defs]
		self.y_names=_y_defs
		self.d_names=[name for name,_,_ in _d_defs]
		self.x_names=[name for name,default,dist,mask in _x_defs]
	
	def eta(self, theta, d, x=[], err=True):
		if x == []: #default x
			x = self.x_default
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition eta: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition eta: d size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#apply discrete mask to theta and d and x
		theta_masked = [(math.floor(tt) if self.theta_masks[i]=='discrete' else tt) for i,tt in enumerate(theta)]
		d_masked = [(math.floor(dd) if self.d_masks[i]=='discrete' else dd) for i,dd in enumerate(d)]
		#make dicts for inputs w/ defs, to ensure consistency
		theta_dict = dict(zip(self.theta_names, theta_masked))
		d_dict = dict(zip(self.d_names, d_masked))
		x_dict = dict(zip(self.x_names, x))
		return self._internal_eta(theta_dict, d_dict, x_dict, err)
	
	def G(self, d, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition G: theta size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition G: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#apply discrete mask to d and x
		d_masked = [(math.floor(dd) if self.d_masks[i]=='discrete' else dd) for i,dd in enumerate(d)]
		#make dicts for inputs w/ defs, to ensure consistency
		d_dict = dict(zip(self.d_names, d_masked))
		x_dict = dict(zip(self.x_names, x))
		return self._internal_G(d_dict, x_dict)
	
	def H(self, theta, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition H: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition H: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#apply discrete mask to theta and x
		theta_masked = [(math.floor(tt) if self.theta_masks[i]=='discrete' else tt) for i,tt in enumerate(theta)]
		#make dicts for inputs w/ defs, to ensure consistency
		theta_dict = dict(zip(self.theta_names, theta_masked))
		x_dict = dict(zip(self.x_names, x))
		return self._internal_H(theta_dict, x_dict)
	
	
	_allowable_prior_types = {
		'gaussian': 2, #mu, sigma
		'gamma_ab': 2, #alpha, beta
		'gamma_mv': 2, #mean, variance
		'lognorm': 2, #mean, variance
		'uniform': 2, #left, right
		'nonrandom': 1, #return value
		'gp_expquad': 4 #variance, ls, prior_pts, mean_fn
		} #knot list (x,y), order, y_stddev, domain, range
		
	
	def prior_rvs(self, num_vals):
		return self._dist_rvs(num_vals, self.priors)
		
	def sample_d(self, num_vals):
		d = self._dist_rvs(num_vals, self.d_dists)
		d_masked = [(math.floor(dd) if self.d_masks[i]=='discrete' else dd) for i,dd in enumerate(d)]
		return d_masked
		
	def sample_x(self, num_vals):
		return self._dist_rvs(num_vals, self.x_dists)
	
	def _dist_rvs(self, num_vals, dist_def):
		vals = [] #a list length num_vals of random numbers of size dim_theta
		for prior in dist_def: ###iterate over dim_theta
			dtype = prior[0]
			params = prior[1]
			thetas_i = [] #need to do this carefully, we have multiple thetas and multiple samples
	
			#generate the rvs for this one particular theta
			if dtype == 'gaussian':
				mu = params[0]
				sigma = params[1]
				thetas_i = scipy.stats.norm.rvs(size=num_vals, loc=mu, scale=sigma)
			elif dtype == 'gamma_ab':
				alpha = params[0]
				beta = params[1]
				thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha, scale=1.0/beta)
			elif dtype == 'gamma_mv':
				mean = params[0]
				variance = params[1]
				alpha = mean**2 / variance
				beta = mean / variance
				thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha, scale=1.0/beta)
			elif dtype == 'lognorm':
				mu = params[0]
				sigma = params[1]
				thetas_i = scipy.stats.lognorm.rvs(size=num_vals, s=sigma, scale=np.exp(mu))
			elif dtype == 'uniform':
				left = params[0]
				right = params[1]
				thetas_i = scipy.stats.uniform.rvs(size=num_vals, loc=left, scale=right-left) #dist is [loc, loc + scale]
			elif dtype == 'nonrandom':
				thetas_i = [param[0] for _ in range(num_vals)]
			elif dtype == 'gp_expquad':
				variance = params[0]
				ls = params[1]
				prior_pts = params[2]
				mean_fn = params[3]
				thetas_i = [sample_gp_prior(variance, ls, prior_pts, mean_fn) for _ in range(num_vals)]
			else:
				raise ValueError("_dist_rvs did not expect prior type "+str(type))
				
			vals.append(thetas_i)
	
		#turn from list of rvs at each prior, to list of theta rvs
		vals = np.transpose(vals)
	
		#return
		if len(vals) == 1:
			return vals[0] #this is a list (dim_theta)
		else:
			return vals #this is a list (num_vals) of random variables (dim_theta)
	
	#assumes single theta, not list of thetas
	"""
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
			elif type == 'gamma_ab':
				alpha = params[0]
				beta = params[1]
				prob_i = scipy.stats.gamma.pdf(theta_i, a=alpha, scale=1.0/beta)
			elif type == 'gamma_mv':
				mean = params[0]
				variance = params[1]
				alpha = mean**2 / variance
				beta = mean / variance
				prob_i = scipy.stats.gamma.pdf(theta_i, a=alpha, scale=1.0/beta)
			elif type == 'lognorm':
				mu = params[0]
				sigma = params[1]
				prob_i = scipy.stats.lognorm.pdf(theta_i, s=sigma, scale=np.exp(mu))
			elif type == 'nonrandom':
				prob_i = (theta_i == params[0])
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
	"""
	
	def theta_nominal(self):
		tnom = []
		for prior in self.priors: ###iterate over dim_theta
			dtype = prior[0]
			params = prior[1]
	
			#generate the rvs for this one particular theta
			if dtype == 'gaussian':
				mu = params[0]
				tnom.append(mu)
			elif dtype == 'gamma_ab':
				alpha = params[0]
				beta = params[1]
				tnom.append(alpha/beta)
			elif dtype == 'gamma_mv':
				mean = params[0]
				tnom.append(mean)
			elif dtype == 'lognorm':
				mu = params[0]
				tnom.append(mu)
			elif dtype == 'uniform':
				left = params[0]
				right = params[1]
				tnom.append(right-left)
			elif dtype == 'nonrandom':
				tnom.append(param[0])
			elif dtype == 'gp_expquad':
				prior_pts = params[2]
				mean_fn = params[3]
				funct = define_functional_mean(prior_pts, mean_fn)
				tnom.append(funct)
			else:
				return 0 #fail!
		return tnom
		
	def __str__(self): #handsome little format for printing the object
		printout = "ProblemDefinition printout:\n"
		printout += self.__repr__() + '\n'
		for key,val in vars(self).items():
			printout += str(key) + " : " + str(val) + '\n'
		return printout


if __name__ == "__main__":
    #unit testing
	def eta(_theta, _d, _x):
		d1 = _d["d1"]
		d2 = _d["d2"]
		t1 = _theta["theta1"]
		t2 = _theta["theta2"]
		factor = _x["y2 factor"]
	
		y = [0,0]
		y[0] =  d1*t1 + t2 + scipy.stats.norm.rvs(loc=0, scale=d2, size=1)[0]
		y[1] = factor*t1/t2 + scipy.stats.norm.rvs(loc=0, scale=d1, size=1)[0]
		return y
		
	def H(_theta, _x):
		t1 = _theta["theta1"]
		t2 = _theta["theta2"]
		K = _x["H factor"]
		
		return K*t1 + t2**2
		
	def Gamma(_d, _x):
		d1 = _d["d1"]
		d2 = _d["d2"]
		C = _x["Cost factor"]
		
		return C/(d1+d2)				   
	
	toy_theta_defs = [ 
						["theta1", ["uniform", [-2,-1]], "continuous"],
						["theta2", ["uniform", [1, 21]], "continuous"]
					 ]
	
	toy_y_defs = ["y1", "y2"]
	
	toy_d_defs = 	[
						["d1", ["uniform", [1, 10]], "discrete" ],
						["d2", ["uniform", [0.01, 1]], "continuous" ]
					]
	
	toy_x_defs = [
					["y2 factor", ["uniform", [.1, 2]], "continuous", 1],
					["H factor", ["uniform", [1, 12]], "discrete", 10],
					["Cost factor", ["uniform", [5, 15]], "continuous", 10]
				 ]

	
	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	toy = ProblemDefinition(eta, H, Gamma, toy_theta_defs, toy_y_defs, toy_d_defs, toy_x_defs)
	print(toy)
	"""print(toy.prior_pdf_unnorm([-1.5,1.5]))"""
	print(toy.prior_rvs(5))
	print(toy.prior_rvs(1))
	print(toy.sample_d(1))