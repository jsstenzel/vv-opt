import sys
import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

import chaospy

sys.path.append('..')
from uq.plotmatrix import *

"""
Taking a similar approach as detailed in Huan & Marzouk 2013
But not using dimension-adaptive sparse quadrature
N: quadrature order
M: polynomial order, the truncation
"""
def likelihood_G_PCE_psp(problem, N, M):
	###Start with problem formulation
	#https://chaospy.readthedocs.io/en/master/reference/distribution/collection.html
	prob_dists = []
	for dist in problem.priors+problem.d_dists:
		if dist[0]=="gaussian": #mu, sigma
			prob_dists.append(chaospy.Normal(dist[1][0],dist[1][1]))
		elif dist[0]=="uniform": #left, right
			prob_dists.append(chaospy.Uniform(dist[1][0],dist[1][1]))
		elif dist[0]=='gamma_ab': #alpha, beta
			k = dist[1][0]
			theta = 1/dist[1][1]
			prob_dists.append(chaospy.Gamma(k,theta))
		elif dist[0]=='gamma_mv': #mean, variance
			m = dist[1][0]
			v = dist[1][1]
			k = m**2 / v
			theta = v / m
			prob_dists.append(chaospy.Gamma(k,theta))
		elif dist[0]=='lognorm': #mean, sigma
			prob_dists.append(chaospy.Lognormal(mu=dist[1][0],sigma=dist[1][1]))
		else:
			print("Unknown distribution for chaospy:",dist)
			exit()
	
	print("Making joint distribution...",flush=True)
	joint = chaospy.J(*prob_dists) #https://chaospy.readthedocs.io/en/master/api/chaospy.J.html?highlight=chaospy.J
	
	def model_solver(parameters):
		#unpack theta,d
		#print(parameters)
		theta = parameters[:problem.dim_theta]
		d = parameters[problem.dim_theta:]
		#print(theta)
		#print(d)
		#exit()
		return problem.G(theta,d)
	
	##this is just to evaluate the joint distribution you defined:
	#parameter_samples = joint.sample(10000)
	#model_evaluations = numpy.array([model_solver(sample) for sample in parameter_samples.T])
	
	###Generate quadrature
	print("Generating quadrature...",flush=True)
	quads = [
		chaospy.generate_quadrature(
		order, #generated for each order up to N
		joint, #pass in the joint distribution we constructed
		sparse=True, #use Smolyak construction
		rule=None) #ensures Clenshaw-Curis quadrature
		for order in range(1, N+1)
	]

	###Evaluate the model
	print("Generating model evals...",flush=True)
	evals = [
		np.array([model_solver(node) for node in nodes.T])
		for nodes, weights in quads
	]
	
	for eval_ in evals:
		print(len(eval_))

	###Psuedo-spectral projection
	print("Generating expansions...",flush=True)
	expansions = [chaospy.generate_expansion(order, joint) for order in range(1, M+1)]
	#expansions[0].round(10)

	#print(len(expansions))
	#print(len(quads))
	#print(len(evals))

	###Generate PCE model
	#IM DOING SOMETHING WRONG HERE - why are we zipping together M expansions and N quadratures?
	print("Generating PCE model...",flush=True)
	model_approx = [
		chaospy.fit_quadrature(expansion, nodes, weights, evals)
		for expansion, (nodes, weights), evals in zip(expansions, quads, evals)
	]
	
	for thing in model_approx:
		print(len(thing))

	expected = chaospy.E(model_approx, joint)
	std = chaospy.Std(model_approx, joint)
	print(expected, flush=True)
	print(std, flush=True)

	return model_approx, joint
	
def likelihood_G_PCE_linregress(problem, n_eval, order):
	###Start with problem formulation
	#https://chaospy.readthedocs.io/en/master/reference/distribution/collection.html
	prob_dists = []
	for dist in problem.priors+problem.d_dists:
		if dist[0]=="gaussian": #mu, sigma
			prob_dists.append(chaospy.Normal(dist[1][0],dist[1][1]))
		elif dist[0]=="uniform": #left, right
			prob_dists.append(chaospy.Uniform(dist[1][0],dist[1][1]))
		elif dist[0]=='gamma_ab': #alpha, beta
			k = dist[1][0]
			theta = 1/dist[1][1]
			prob_dists.append(chaospy.Gamma(k,theta))
		elif dist[0]=='gamma_mv': #mean, variance
			m = dist[1][0]
			v = dist[1][1]
			k = m**2 / v
			theta = v / m
			prob_dists.append(chaospy.Gamma(k,theta))
		elif dist[0]=='lognorm': #mean, sigma
			prob_dists.append(chaospy.Lognormal(mu=dist[1][0],sigma=dist[1][1]))
		else:
			print("Unknown distribution for chaospy:",dist)
			exit()
	
	print("Making joint distribution...",flush=True)
	joint = chaospy.J(*prob_dists) #https://chaospy.readthedocs.io/en/master/api/chaospy.J.html?highlight=chaospy.J

	print("Sampling from d,theta space...", flush=True)
	samples = joint.sample(n_eval, rule="sobol")
	
	print("Evaluating...", flush=True)
	G_evals = np.array([problem.G(sample[:problem.dim_theta],sample[problem.dim_theta:]) for sample in samples.T])

	print("Fitting...", flush=True)
	expansion = chaospy.generate_expansion(order=order, dist=joint)
	proxy_model, coeffs = chaospy.fit_regression(expansion, samples, G_evals, retall=True)
	
	return proxy_model

	

	
	
def sample_likelihood_pce(problem, gauss_model_approx):
	stochastic_dim = problem.dim_theta + problem.dim_d
	total_multindex = totalOrderMultiIndices(stochastic_dim, p)

def eval_likelihood_pce(problem, gauss_model_approx):
	stochastic_dim = problem.dim_theta + problem.dim_d
	total_multindex = totalOrderMultiIndices(stochastic_dim, p)
	


	
def error_in_mean(predicted_mean, true_mean):
	"""
	How close the estimated mean is the to the true mean.

	Args:
		predicted_mean (numpy.ndarray):
			The estimated mean.
		true_mean (numpy.ndarray):
			The reference mean value. Must be same shape as
			``prediction_mean``.

	Returns:
		(float):
			The mean absolute distance between predicted
			and true values.
	"""
	return numpy.mean(numpy.abs(predicted_mean-true_mean))

def error_in_variance(predicted_variance,
					  true_variance):
	"""
	How close the estimated variance is the to the true variance.

	Args:
		predicted_variance (numpy.ndarray):
			The estimated variance.
		true_variance (numpy.ndarray):
			The reference variance value.
			Must be same shape as ``predicted_variance``.

	Returns:
		(float):
			The mean absolute distance between
			predicted and true values.
	"""
	return numpy.mean(numpy.abs(predicted_variance-true_variance))