#pdf estimation methods, including kernel density estimations

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

sys.path.append('../..')
#from approx.regression_models import *
from problems.problem_definition import *
from uq.plotmatrix import *
from inference.goal_based_inference import *

#train a model on the data p(theta),d -> y and p(theta) -> Q
#to develop the joint model Q = J(y,d)
#for use in the goal-based approach
#with simple OLS linear regression!
def joint_model_linear(problem, N, doPrint=False, doPlot=False):
	if doPrint:
			print("Generating the training data...", flush=True)
			
	#model inputs
	theta_train = problem.prior_rvs(N)
	d_train = [d for d in problem.sample_d(N)]
	
	#model outputs
	qoi_train = [problem.H(theta) for theta in theta_train]
	y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
	
	#print(theta_train)
	#print(d_train)
	#print(qoi_train)
	#print(y_train)
	
	joint_input = [yi+di for yi,di in zip(y_train, d_train)]
	joint_output = np.array(qoi_train)
	
	#print(joint_input)
	#print(joint_output)
	
	if doPrint:
			print("Training the joint model...", flush=True)

	# Create linear regression object
	regr = linear_model.LinearRegression()

	# Train the model using the training sets
	regr.fit(joint_input, joint_output)

	# Make self-predictions using the training set, to compare to actual outputs
	joint_pred = regr.predict(joint_input)

	# Print the joint model parameters
	print("Coefficients: \n", regr.coef_)
	print("MSE of training data: %.2f" % mean_squared_error(joint_pred, joint_output))
	# The coefficient of determination: 1 is perfect prediction
	print("Coefficient of determination: %.2f" % r2_score(joint_pred, joint_output))

	# Plot outputs?
	#if doPlot:
	#	plotfn_large_domain()
	
	model = regr.predict
	MSE = mean_squared_error(joint_pred, joint_output)
	
	return model, MSE
	
#train a model on the data p(theta),d -> y and p(theta) -> Q
#to develop the joint model Q = J(y,d)
#for use in the goal-based approach
#with the same GMM approach I used originally
def joint_model_gmm(problem, N, doPrint=False, doPlot=False):
	if doPrint:
			print("Generating the training data...", flush=True)
			
	#model inputs
	theta_train = problem.prior_rvs(N)
	d_train = [d for d in problem.sample_d(N)]
	
	#model outputs
	qoi_train = [problem.H(theta) for theta in theta_train]
	y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
	
	#print(theta_train)
	#print(d_train)
	#print(qoi_train)
	#print(y_train)
	
	joint_input = [yi+di for yi,di in zip(y_train, d_train)]
	joint_output = np.array(qoi_train)
	
	#print(joint_input)
	#print(joint_output)
	
	if doPrint:
			print("Training the joint model...", flush=True)

	# Create linear regression object
	gmm = gbi_train_model(joint_output, joint_output, joint_input, ncomp=0, verbose=doPrint)
	
	model = gmm
	
	return model

def polynomial_features(X, k):
	return np.hstack([X ** i for i in range(1, k + 1)])

#train a model on the data p(theta),d -> y and p(theta) -> Q
#to develop the joint model Q = J(y,d)
#for use in the goal-based approach
#using Bayesian Regression with Polynomial Features
def joint_model_bayes_polynomial(problem, N, pk, doPrint=False, doPlot=False):
	0
	#TODO
	
	
def joint_model_PCE_linregress(problem, N, order, sobol=False):
	if doPrint:
			print("Generating the training data...", flush=True)
			
	if sobol == False:
		#model inputs
		theta_train = problem.prior_rvs(N)
		d_train = [d for d in problem.sample_d(N)]
		
		#model outputs
		qoi_train = [problem.H(theta) for theta in theta_train]
		y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
		
		#print(theta_train)
		#print(d_train)
		#print(qoi_train)
		#print(y_train)
		
		joint_input = [yi+di for yi,di in zip(y_train, d_train)]
		joint_output = np.array(qoi_train)
		
		#print(joint_input)
		#print(joint_output)
		
	else:
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
		
		if doPrint:
			print("Making joint distribution...",flush=True)
		joint = chaospy.J(*prob_dists) #https://chaospy.readthedocs.io/en/master/api/chaospy.J.html?highlight=chaospy.J

		if doPrint:
			print("Sampling from d,theta space...", flush=True)
		samples = joint.sample(n_eval, rule="sobol")
		
		#separate the inputs
		theta_train = [sample[:problem.dim_theta] for sample in samples.T]
		d_train = [sample[problem.dim_theta:] for sample in samples.T]
		
		#model outputs
		qoi_train = [problem.H(theta) for theta in theta_train]
		y_train = [problem.eta(theta, d) for theta,d in zip(theta_train,d_train)]
		
		joint_input = [yi+di for yi,di in zip(y_train, d_train)]
		joint_output = np.array(qoi_train)
	
	if doPrint:
			print("Training the joint model...", flush=True)
			
	###no no, this is wrong, I need to make a joint distribution

	expansion = chaospy.generate_expansion(order=order, dist=joint)
	proxy_model, coeffs = chaospy.fit_regression(expansion, joint_input, joint_output, retall=True)
	
	return proxy_model