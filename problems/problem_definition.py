import sys
import os
import scipy.stats
import scipy.interpolate
import numpy as np
import math
from copy import deepcopy

sys.path.append('..')
from approx.gaussian_process import *


class ProblemDefinition:
	def __init__(self, _eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs):
		self._internal_eta = _eta #eta(theta, d, x)
		self._internal_H = _H #H(theta, x)
		self._internal_G = _G #G(d, x)
		
		#First, need to process the multivariate priors in theta
		thetanames=[]
		thetamasks=[]
		for theta_def in _theta_defs:
			name = theta_def[0]
			dtype = theta_def[1][0]
			params = theta_def[1][1]
			mask = theta_def[2]
			if dtype == "gaussian_multivar":
				thetanames.extend([name+str(i) for i,_ in enumerate(params[0])])
				thetamasks.extend([mask for _ in params[0]])
			else:
				thetanames.append(name)
				thetamasks.append(mask)
		dimtheta=len(thetanames)
		
		#get dimensions, set default parameters
		self.dim_d = len(_d_defs)
		self.dim_y = len(_y_defs)
		self.dim_x = len(_x_defs)
		self.dim_theta = dimtheta #self.dim_theta = len(_theta_defs)
		self.x_default = [default for name,dist,mask,default in _x_defs]
		
		#hack the x_defs to save wasted effort:
		for xdef in _x_defs:
			if xdef[1] == []:
				xdef[1] = ["nonrandom", [xdef[3]]]
		
		#priors is a list of defintions of prior distributions on the vector of thetas
		#so it's a list of pairs, first is type and second is the list of corresponding parameters
		#type checks:
		for _,prior,mask in _theta_defs + _d_defs + [x[:-1] for x in _x_defs]:
			dtype = prior[0]
			params = prior[1]
			
			if dtype not in self._allowable_prior_types.keys():
				raise ValueError("Incorrect prior probability function definition: "+str(dtype)+" not recognized.")
			if len(params) != self._allowable_prior_types[dtype]:
				raise ValueError('Wrong number of arguments for prior type '+str(dtype)+'; got '+str(len(params))+' and expected '+str(self._allowable_prior_types[dtype]))
			if not (mask in ['continuous','discrete','functional','object','discretized100']):
				raise ValueError('Variable mask not an expected value: '+str(mask))

		#if you pass all of that,
		self.priors  = [prior for _,prior,_ in _theta_defs] #now this is <= dim_theta w/ gaussian_multivar
		self.d_dists = [dist for _,dist,_ in _d_defs]
		self.x_dists = [dist for _,dist,_,_ in _x_defs]
		
		self.theta_masks = thetamasks #self.theta_masks = [mask for _,_,mask in _theta_defs]
		self.d_masks     = [mask for _,_,mask in _d_defs]
		self.x_masks     = [mask for _,_,mask,_ in _x_defs]
		
		#Get theta_nominal only once, to save repetition
		self.theta_nominal = self._theta_nominal()
		
		#for documentation:
		self.theta_names= thetanames #self.theta_names=[name for name,_,_ in _theta_defs]
		self.y_names=_y_defs
		self.d_names=[name for name,_,_ in _d_defs]
		self.x_names=[name for name,default,dist,mask in _x_defs]
		
		#pre-construct some of the dicts, so I never have to repeat it unnecessarily
		self.x_dict = dict(zip(self.x_names, self.x_default))
		self.prior_mean_dict = dict(zip(self.theta_names, self.theta_nominal))
	
	def eta(self, theta, d, x=[], err=True):
		if x == []: #default x
			x_dict = deepcopy(self.x_dict)
		else:
			x_dict = dict(zip(self.x_names, x))
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition eta: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition eta: d size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x_dict) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#apply discrete mask to theta and d and x
		theta_masked = self._mask(theta, self.theta_masks, self.priors)
		d_masked = self._mask(d, self.d_masks, self.d_dists)
		#make dicts for inputs w/ defs, to ensure consistency
		theta_dict = dict(zip(self.theta_names, theta_masked))
		d_dict = dict(zip(self.d_names, d_masked))
		#provide prior means here, in case they're needed for imputation
		return self._internal_eta(theta_dict, d_dict, x_dict, self.prior_mean_dict, err)
	
	def G(self, d, x=[]):
		if x == []: #default x
			x_dict = deepcopy(self.x_dict)
		else:
			x_dict = dict(zip(self.x_names, x))
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition G: theta size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x_dict) != self.dim_x:
			raise ValueError("Input to ProblemDefinition G: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#apply discrete mask to d and x
		d_masked = self._mask(d, self.d_masks, self.d_dists)
		#make dicts for inputs w/ defs, to ensure consistency
		d_dict = dict(zip(self.d_names, d_masked))
		return self._internal_G(d_dict, x_dict)
	
	def H(self, theta, x=[], verbose=False):
		if x == []: #default x
			x_dict = deepcopy(self.x_dict)
		else:
			x_dict = dict(zip(self.x_names, x))
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition H: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(x_dict) != self.dim_x:
			raise ValueError("Input to ProblemDefinition H: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#apply discrete mask to theta and x
		theta_masked = self._mask(theta, self.theta_masks, self.priors)
		#make dicts for inputs w/ defs, to ensure consistency
		theta_dict = dict(zip(self.theta_names, theta_masked))
		return self._internal_H(theta_dict, x_dict, verbose)
		
	def _mask(self, vec, masks, defs):
		masked_vec = [None]*len(vec)
		for i,dd in enumerate(vec):
			if masks[i]=='discrete':
				masked_vec[i] = math.floor(dd)
			elif masks[i]=='discretized100':
				#discretizes the min-max space to 100 equal points, then rounds dd to nearest point
				ddomain = np.linspace(float(defs[i][1][0]), float(defs[i][1][1]), 101)
				masked_vec[i] = min(ddomain, key=lambda y: abs(dd - y))
			else:
				masked_vec[i] = dd
		return masked_vec
	
	_allowable_prior_types = {
		'gaussian': 2, #mu, sigma
		'gaussian_multivar': 2, #mean vector, covariance
		'gamma_ab': 2, #alpha, beta
		'gamma_mv': 2, #mean, variance
		'beta': 2, #a, b
		'lognorm': 2, #mean, variance
		'uniform': 2, #left, right
		'nonrandom': 1, #return value
		'gp_expquad': 4 #variance, ls, prior_pts, mean_fn
		} #knot list (x,y), order, y_stddev, domain, range
		
	
	def prior_rvs(self, num_vals):
		return self._dist_rvs(num_vals, self.priors)
		
	def sample_d(self, num_vals):
		d = self._dist_rvs(num_vals, self.d_dists)
		if num_vals == 1: #stupid, inefficient
			d = [d]
		d_masked = [self._mask(dj, self.d_masks, self.d_dists) for dj in d]
		
		if num_vals == 1: #also somewhat stupid
			return d_masked[0] #this is a list (dim_d)
		else:
			return d_masked #this is a list (num_vals) of random variables (dim_d)
		
	def sample_x(self, num_vals):
		return self._dist_rvs(num_vals, self.x_dists)
	
	def _dist_rvs(self, num_vals, dist_def):
		vals = [] #a list length num_vals of random numbers of size dim_theta
		for prior in dist_def: ###iterate over dim_theta
			dtype = prior[0]
			params = prior[1]
			#need to do this carefully, we have multiple thetas and multiple samples
	
			#generate the rvs for this one particular theta
			if dtype == 'gaussian':
				mu = params[0]
				sigma = params[1]
				thetas_i = scipy.stats.norm.rvs(size=num_vals, loc=mu, scale=sigma)
				vals.append(thetas_i.tolist())
			elif dtype == 'gaussian_multivar':
				mean_vector = np.array(params[0])
				covariance = np.array(params[1])
				multisamples = scipy.stats.multivariate_normal.rvs(size=num_vals, mean=mean_vector, cov=covariance)
				for thetas_i in multisamples.T: #i think im breaking this down in the wrong direction...
					vals.append(thetas_i.tolist())
			elif dtype == 'gamma_ab':
				alpha = params[0]
				beta = params[1]
				thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha, scale=1.0/beta)
				vals.append(thetas_i.tolist())
			elif dtype == 'gamma_mv':
				mean = params[0]
				variance = params[1]
				alpha = mean**2 / variance
				beta = mean / variance
				thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha, scale=1.0/beta)
				vals.append(thetas_i.tolist())
			elif dtype == 'beta':
				a = params[0]
				b = params[1]
				thetas_i = scipy.stats.beta.rvs(a=a, b=b, size=num_vals)
				vals.append(thetas_i.tolist())
			elif dtype == 'lognorm':
				mu = params[0]
				sigma = params[1]
				thetas_i = scipy.stats.lognorm.rvs(size=num_vals, s=sigma, scale=np.exp(mu))
				vals.append(thetas_i.tolist())
			elif dtype == 'uniform':
				left = params[0]
				right = params[1]
				thetas_i = scipy.stats.uniform.rvs(size=num_vals, loc=left, scale=right-left) #dist is [loc, loc + scale]
				vals.append(thetas_i.tolist())
			elif dtype == 'nonrandom':
				thetas_i = [params[0] for _ in range(num_vals)]
				vals.append(thetas_i)
			elif dtype == 'gp_expquad':
				variance = params[0]
				ls = params[1]
				prior_pts = params[2]
				mean_fn = params[3]
				#TODO to make theta-sampling faster, I should save these in the problem
				gp_prior = GaussianProcessDist1D(variance, ls, prior_pts, mean_fn)
				thetas_i = [gp_prior.sample() for _ in range(num_vals)]
				vals.append(thetas_i)
			else:
				raise ValueError("_dist_rvs did not expect prior type "+str(type))

		#turn from list of rvs at each prior, to list of theta rvs
		vals = np.transpose(vals)
	
		#return
		if len(vals) == 1:
			return vals[0] #this is a list (dim_theta)
		else:
			return vals #this is a list (num_vals) of random variables (dim_theta)
			
	def print_theta(self, theta):
		for i,prior in enumerate(self.priors):
			dtype = prior[0]
			params = prior[1]
			print(self.theta_names[i],theta[i], flush=True)
			if dtype == 'gp_expquad':
				gp_sample = theta[i]
				gp_sample.plot_prior()
	
	#assumes single theta, not list of thetas

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
	
	def _theta_nominal(self):
		"""
		Standard behavior is to return a list of nominal error-free values of every theta, including error-free mean functions
		This is so you can get a "nominal" value for y with the following expression:
		y_nominal = problem.eta(problem.theta_nominal(), d, err=False)
		
		For GP's, the behavior is to provide a GP that is exactly evaluated at the value of the mean function
		"""
		tnom = []
		for prior in self.priors: ###iterate over dim_theta
			dtype = prior[0]
			params = prior[1]
	
			#generate the rvs for this one particular theta
			if dtype == 'gaussian':
				mu = params[0]
				tnom.append(mu)
			elif dtype == 'gaussian_multivar':
				mean_vector = params[0]
				for mu in mean_vector:
					tnom.append(mu)
			elif dtype == 'gamma_ab':
				alpha = params[0]
				beta = params[1]
				tnom.append(alpha/beta)
			elif dtype == 'gamma_mv':
				mean = params[0]
				tnom.append(mean)
			elif dtype == 'beta':
				a = params[0]
				b = params[1]
				tnom.append(a/(a+b))
			elif dtype == 'lognorm':
				mu = params[0]
				tnom.append(mu)
			elif dtype == 'uniform':
				left = params[0]
				right = params[1]
				tnom.append((right+left)/2.0)
			elif dtype == 'nonrandom':
				tnom.append(param[0])
			elif dtype == 'gp_expquad':
				prior_pts = params[2]
				mean_fn = params[3]
				funct = define_functional_from_meanfn(prior_pts, mean_fn)
				tnom.append(funct)
			else:
				raise ValueError("_theta_nominal couldn't construct vector with "+str(dtype))
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
						["d2", ["uniform", [0.01, 1]], "continuous" ],
						["d3", ["uniform", [10, 20]], "discretized100" ]
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
	print(toy.sample_d(3))
