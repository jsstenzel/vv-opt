import sys
import os
import scipy.stats
import scipy.interpolate
import numpy as np

sys.path.append('..')
from problems.functionals import *


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
		self.x_default = [default for name,default in _x_defs]
		
		#priors is a list of defintions of prior distributions on the vector of thetas
		#so it's a list of pairs, first is type and second is the list of corresponding parameters
		#type checks:
		for _,prior in _theta_defs:
			type = prior[0]
			params = prior[1]
			
			if type not in self._allowable_prior_types.keys():
				raise ValueError("Incorrect prior probability function definition: "+str(type)+" not recognized.")
			if len(params) != self._allowable_prior_types[type]:
				raise ValueError('Wrong number of arguments for prior type '+str(type)+'; got '+str(len(params))+' and expected '+str(self._allowable_prior_types[type]))
			#for param in params:
			#	if not isinstance(param, (int, float)):
			#		raise ValueError('Wrong kind of arguments for prior type '+str(type)+', need int or float.')

		#if you pass all of that,
		self.priors = [prior for _,prior in _theta_defs]
	
		#for documentation:
		self.theta_names=[name for name,_ in _theta_defs]
		self.y_names=_y_defs
		self.d_names=_d_defs
		self.x_names=[name for name,_ in _x_defs]
	
	def eta(self, theta, d, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition eta: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition eta: d size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#make dicts for inputs w/ defs, to ensure consistency
		theta_dict = dict(zip(self.theta_names, theta))
		d_dict = dict(zip(self.d_names, d))
		x_dict = dict(zip(self.x_names, x))
		return self._internal_eta(theta_dict, d_dict, x_dict)
	
	def G(self, theta, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(theta) != self.dim_theta:
			raise ValueError("Input to ProblemDefinition eta: theta size "+str(self.dim_theta)+" expected, "+str(len(theta))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#make dicts for inputs w/ defs, to ensure consistency
		theta_dict = dict(zip(self.theta_names, theta))
		x_dict = dict(zip(self.x_names, x))
		return self._internal_G(theta_dict, x_dict)
	
	def H(self, d, x=[]):
		if x == []: #default x
			x = self.x_default
		if len(d) != self.dim_d:
			raise ValueError("Input to ProblemDefinition eta: d size "+str(self.dim_d)+" expected, "+str(len(d))+" provided.")
		if len(x) != self.dim_x:
			raise ValueError("Input to ProblemDefinition eta: x size "+str(self.dim_x)+" expected, "+str(len(x))+" provided.")
		#make dicts for inputs w/ defs, to ensure consistency
		d_dict = dict(zip(self.d_names, d))
		x_dict = dict(zip(self.x_names, x))
		return self._internal_H(d_dict, x_dict)
	
	
	_allowable_prior_types = {
		'gaussian': 2, #mu, sigma
		#'pos_gaussian': 2, #mu, sigma
		#'bound_gaussian': 4 #mu, sigma, left, right
		'uniform': 2,
		'funct_splines': 5} #knot list (x,y), order, y_stddev, domain, range
		
	
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
			#if type == 'pos_gaussian':
			#	mu = params[0]
			#	sigma = params[1]
			#	thetas_i=0
			#	while thetas_i <= 0:
			#		thetas_i = scipy.stats.norm.rvs(size=num_vals, loc=mu, scale=sigma) #fix num_vals thing
			#if type == 'bound_gaussian':
			#	mu = params[0]
			#	sigma = params[1]
			#	thetas_i = scipy.stats.norm.rvs(size=num_vals, loc=mu, scale=sigma)
			#	while thetas_i < params[2] or thetas_i > params[3]:
			#		thetas_i = scipy.stats.norm.rvs(size=num_vals, loc=mu, scale=sigma) #fix num_vals thing
			elif type == 'uniform':
				left = params[0]
				right = params[1]
				thetas_i = scipy.stats.uniform.rvs(size=num_vals, loc=left, scale=right-left) #dist is [loc, loc + scale]
			elif type == 'funct_splines':
				for _ in range(num_vals):
					data = params[0] #list of pairs of points to define the prior mean
					order = params[1] #k order of spline interpolated fit
					y_stddev = params[2]
					xmin = params[3][0]
					xmax = params[3][1] #domain=(xmin,xmax) boundary on the priors
					ymin = params[4][0]
					ymax = params[4][1] #range=(ymin,ymax) boundary on the priors
					
					x = [point[0] for point in data]
					y = [point[1] + scipy.stats.norm.rvs(size=1, loc=0, scale=y_stddev)[0] for point in data]
					
					sample = Functional(x, y)
					sample.set_xlim(xmin, xmax)
					sample.set_ylim(ymin, ymax)
					sample.spline_interp(order)
					
					thetas_i.append(sample)
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
			#if type == 'pos_gaussian':
			#	mu = params[0]
			#	sigma = params[1]
			#	prob_i = scipy.stats.norm.pdf(theta_i, loc=mu, scale=sigma)
			#	if theta_i <= 0:
			#		prob_i = 0
			#if type == 'bound_gaussian':
			#	mu = params[0]
			#	sigma = params[1]
			#	prob_i = scipy.stats.norm.pdf(theta_i, loc=mu, scale=sigma)
			#	if theta_i <= params[2] or theta_i >= params[3]:
			#		prob_i = 0
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
						("theta1", ["uniform", [-2,-1]]),
						("theta2", ["uniform", [1, 21]])
					 ]
	
	toy_y_defs = ["y1", "y2"]
	
	toy_d_defs = ["d1", "d2"]
	
	toy_x_defs = [("y2 factor", 1),("H factor", 10),("Cost factor", 10)]

	
	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	toy = ProblemDefinition(eta, H, Gamma, toy_theta_defs, toy_y_defs, toy_d_defs, toy_x_defs)
	print(toy)
	print(toy.prior_pdf_unnorm([-1.5,1.5]))
	print(toy.prior_rvs(5))
	print(toy.prior_rvs(1))
