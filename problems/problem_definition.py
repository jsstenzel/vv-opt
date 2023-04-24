import sys
import os
import scipy.stats

sys.path.append('../..')
#from obed.obed import *

class ProblemDefinition:
	def __init__(self, _eta, _H, _G, _prior_type, _prior_params):
		self.eta = _eta
		self.H = _H
		self.G = _G
		#dimension of y?
		#dimension of theta?
		#dimension of d?
		#dimension of x?
		
		if _prior_type in self._allowable_prior_types.keys():
			self.prior_type = _prior_type
		else:
			raise ValueError("Incorrect prior probability function definition:"+str(_prior_type)+" not recognized.")
		
		if len(_prior_params) == self._allowable_prior_types[_prior_type]:
			self.prior_params = _prior_params
		else:
			raise ValueError('Wrong number of arguments for prior type '+str(_prior_type)+'; got '+str(len(_prior_params))+' and expected '+str(self._allowable_prior_types[_prior_type]))
		
	_allowable_prior_types = {'gaussian': 2, #mu, sigma
								   'uniform': 2} #left bound, right bound
	
	#oh these need to be multidimensional. fix that next
	def prior_rvs(self, num_vals):
		#generate
		if self.prior_type == 'gaussian':
			mu = self.prior_params[0]
			sigma = self.prior_params[1]
			vals = scipy.stats.norm.rvs(size=num_vals, loc=mu, scale=sigma)
		elif self.prior_type == 'uniform':
			left = self.prior_params[0]
			right = self.prior_params[1]
			vals = scipy.stats.uniform.rvs(size=num_vals, loc=left, scale=right-left) #dist is [loc, loc + scale]
		else:
			raise ValueError("prior_rvs did not expect prior type "+str(self.prior_type))
			
		#return
		if len(vals) == 1:
			return vals[0]
		else:
			return vals
	
	def prior_pdf(self, theta):
		#evaluate and return theta 
		if self.prior_type == 'gaussian':
			mu = self.prior_params[0]
			sigma = self.prior_params[1]
			return scipy.stats.norm.pdf(theta, loc=mu, scale=sigma)
		elif self.prior_type == 'uniform':
			left = self.prior_params[0]
			right = self.prior_params[1]
			return scipy.stats.uniform.pdf(theta, loc=left, scale=right-left)
		else:
			raise ValueError("prior_pdf did not expect prior type "+str(self.prior_type))