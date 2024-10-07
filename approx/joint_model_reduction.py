#pdf estimation methods, including kernel density estimations

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

sys.path.append('../..')
from approx.regression_models import *
from problems.problem_definition import *
from uq.plotmatrix import *

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
	
	if doPrint:
			print("Training the joint model...", flush=True)
	
	#print(theta_train)
	#print(d_train)
	#print(qoi_train)
	#print(y_train)
	
	joint_input = [yi+di for yi,di in zip(y_train, d_train)]
	joint_output = np.array(qoi_train)
	
	#print(joint_input)
	#print(joint_output)

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