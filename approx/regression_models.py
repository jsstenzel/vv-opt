#pdf estimation methods, including kernel density estimations

import sys
import numpy as np
import sklearn
import sklearn.linear_model
import sklearn.neighbors
import sklearn.neural_network
#import data_gen #need to replace that!
#import utils #need to replace that!
import pandas as pd
from collections import namedtuple
from scipy.optimize import curve_fit

sys.path.append('..')
from problems.problem_definition import *
from uq.plotmatrix import *

## Linear Regression
"""
h_lin = sklearn.linear_model.LinearRegression().fit(X, Y)
utils.plot_data_and_fit(X, Y, h_lin.predict, -20, 20) # Change the last two inputs to plot a different range of x.
"""
#def lin_reg(X_in, y):
#	fit = sklearn.linear_model.LinearRegression().fit(X, Y)
#	return fit

## Nearest-Neighbor Regression
"""
h_knn = sklearn.neighbors.KNeighborsRegressor(n_neighbors=7).fit(X, Y)
utils.plot_data_and_fit(X, Y, h_knn.predict, -20, 20)
"""

## Neural Network
"""
neural_net = sklearn.neural_network.MLPRegressor(hidden_layer_sizes=(50, 50, 50, 50),
                                                 verbose=False,
                                                 max_iter=1000)
h_net = neural_net.fit(X, Y.ravel())
utils.plot_data_and_fit(X, Y, h_net.predict, -11, 20)
"""

## Linear Regression with Polynomial Features
"""
pk = 9
dpf = utils.polynomial_features(X, pk)
h_lin_pf = sklearn.linear_model.LinearRegression().fit(dpf, Y)
utils.plot_data_and_fit(X, Y,
                        lambda x: h_lin_pf.predict(utils.polynomial_features(x, pk)), -10, 20)
"""

## Linear regression with fourier features
def fourier_features(X, k, div):
    return np.hstack([np.cos(np.pi * i * X / div) for i in range(1, k + 1)] +
                     [np.sin(np.pi * i * X / div) for i in range(1, k + 1)])

"""
ffk = 5
div = 10 # roughly the domain of the data, to ensure we have a low-frequency component
ff = utils.fourier_features(X, ffk, div)
h_lin_ff = sklearn.linear_model.LinearRegression().fit(ff, Y)
utils.plot_data_and_fit(X, Y,
                        lambda x: h_lin_ff.predict(utils.fourier_features(x, ffk, div)), -20, 20)
"""
def linreg_fourier_throughput_file(thru_file, order, bandpass, doPlot=False, doErr=False):
	lambda_pts = []
	thru_pts = []
	with open(thru_file, "r") as f:
		for line in f:
			words = line.split()
			if len(words)==0:
				continue
			if "#" not in words[0]:
				lambda_pts.append(float(words[0]))
				thru_pts.append(float(words[1]))
				
	return linreg_fourier_throughput(lambda_pts, thru_pts, order, bandpass, doPlot=doPlot, doErr=doErr)

def linreg_fourier_throughput(lambda_pts, thru_pts, order, bandpass, doPlot=False, doErr=False):
	domain_min = min(lambda_pts)
	domain_max = max(lambda_pts)

	X=np.array([[pt] for pt in lambda_pts])
	Y=np.array(thru_pts)
	#print(X)
	div = bandpass # roughly the domain of the data, to ensure we have a low-frequency component
	ff = fourier_features(X, order, div)
	#print(ff)
	h_lin_ff = sklearn.linear_model.LinearRegression().fit(ff, Y)
	print(h_lin_ff.coef_, h_lin_ff.intercept_, flush=True)
	
	if doPlot:
		plot_pts = np.array([[pt] for pt in np.linspace(domain_min, domain_max, 100)])
		Yfit = h_lin_ff.predict(fourier_features(plot_pts, order, div))
	
		plt.plot(X, Y, c='k')
		plt.plot(plot_pts, Yfit, c='orange')
		plt.xlim(domain_min, domain_max)
		plt.show()
		
	if doErr:
		Yfit = h_lin_ff.predict(fourier_features(X, order, div))
		Y_diffs = (Yfit - Y)**2
		MSE = np.mean(Y_diffs)
		print("linreg MSE:", MSE)
	
	#print(h_lin_ff.intercept_, h_lin_ff.coef_)
	return h_lin_ff.coef_, h_lin_ff.intercept_, order, div
	
def throughput_from_linfourier_coeffs(coeffs, intercept, order, div, lambda_pts, doPlot=False):
	h_lin_ff = sklearn.linear_model.LinearRegression()
	h_lin_ff.coef_ = np.array(coeffs)
	h_lin_ff.intercept_ = float(intercept)
	
	X=np.array([[pt] for pt in lambda_pts])
	thru_pts = h_lin_ff.predict(fourier_features(X, order, div))
	
	if doPlot:
		domain_min = min(lambda_pts)
		domain_max = max(lambda_pts)
		plot_pts = np.array([[pt] for pt in np.linspace(domain_min, domain_max, 100)])
		Yfit = h_lin_ff.predict(fourier_features(plot_pts, order, div))
	
		plt.plot(plot_pts, Yfit, c='orange')
		plt.scatter(X, thru_pts, c='blue')
		plt.xlim(domain_min, domain_max)
		plt.show()
	
	return thru_pts
	
####################################################

def sigmoid_step(x, lval, step_pt, rval, power):
	scale = abs(rval-lval)
	intercept = min([lval,rval])
	x0 = step_pt
	b = power if lval<=rval else -power

	return scale * scipy.special.expit((x-x0)*b) + intercept
	
def sigmoid_fit_throughput_file(thru_file, doPlot=False, doErr=False):
	lambda_pts = []
	thru_pts = []
	with open(thru_file, "r") as f:
		for line in f:
			words = line.split()
			if len(words)==0:
				continue
			if "#" not in words[0]:
				lambda_pts.append(float(words[0]))
				thru_pts.append(float(words[1]))
				
	return sigmoid_fit_throughput(lambda_pts, thru_pts, doPlot=doPlot, doErr=doErr)
	
def sigmoid_fit_throughput(lambda_pts, thru_pts, doPlot=False, doErr=False):
	domain_min = min(lambda_pts)
	domain_max = max(lambda_pts)
	range_min = min(thru_pts)
	range_max = max(thru_pts)
	
	X=np.array(lambda_pts)
	Y=np.array(thru_pts)
	
	lbound = [range_min, domain_min, range_min, 1/(10*(domain_max-domain_min))]#lower bound of lval, step_pt, rval, power
	rbound = [range_max, domain_max, range_max, (domain_max-domain_min)/len(lambda_pts)]#upper bound of lval, step_pt, rval, power
	popt,pcov=curve_fit(sigmoid_step,X,Y,bounds=(lbound,rbound))
	
	lval = popt[0]
	step_pt = popt[1]
	rval = popt[2]
	power = popt[3]
	
	if doPlot:
		plot_pts = np.array(np.linspace(domain_min, domain_max, len(lambda_pts)*10))
		Yfit = [sigmoid_step(x, lval, step_pt, rval, power) for x in plot_pts]
	
		plt.plot(X, Y, c='k')
		plt.plot(plot_pts, Yfit, c='orange')
		plt.xlim(domain_min, domain_max)
		plt.show()
		
	if doErr:
		Yfit = [sigmoid_step(x, lval, step_pt, rval, power) for x in X]
		Y_diffs = (Yfit - Y)**2
		MSE = np.mean(Y_diffs)
		print("sigmoid MSE:", MSE)
		
	return lval, step_pt, rval, power
	
def throughput_from_sigmoidfit_coeffs(lval, step_pt, rval, power, lambda_pts):
	X=np.array(lambda_pts)
	thru_pts = [sigmoid_step(x, lval, step_pt, rval, power) for x in X]
	
	#plt.plot(X, thru_pts, c='blue')
	#plt.show()
		
	return thru_pts
	
####################################################

def poly_fn(x, params):
	total = 0
	for i,p in enumerate(params):
		total += p * x**(len(params)-1-i)
	return total
	
def poly_fit_throughput_file(thru_file, power, doPlot=False, doErr=False):
	lambda_pts = []
	thru_pts = []
	with open(thru_file, "r") as f:
		for line in f:
			words = line.split()
			if len(words)==0:
				continue
			if "#" not in words[0]:
				lambda_pts.append(float(words[0]))
				thru_pts.append(float(words[1]))
	
	return poly_fit_throughput(lambda_pts, thru_pts, power, doPlot=doPlot, doErr=doErr)
	
def poly_fit_throughput(lambda_pts, thru_pts, power, doPlot=False, doErr=False):
	domain_min = min(lambda_pts)
	domain_max = max(lambda_pts)
	range_min = min(thru_pts)
	range_max = max(thru_pts)
	
	X=np.array(lambda_pts)
	Y=np.array(thru_pts)
	
	popt=np.polyfit(X, Y, power)
	
	if doPlot:
		plot_pts = np.array(np.linspace(domain_min, domain_max, len(lambda_pts)*10))
		Yfit = [poly_fn(x, popt) for x in plot_pts]
	
		plt.plot(X, Y, c='k')
		plt.plot(plot_pts, Yfit, c='orange')
		plt.xlim(domain_min, domain_max)
		plt.show()
		
	if doErr:
		Yfit = [poly_fn(x, popt) for x in X]
		Y_diffs = (Yfit - Y)**2
		MSE = np.mean(Y_diffs)
		print("polyfit MSE:", MSE)
		
	return popt
	
def throughput_from_polyfit_coeffs(popt, lambda_pts):
	X=np.array(lambda_pts)
	thru_pts = [poly_fn(x, popt) for x in X]
	
	#plt.plot(X, thru_pts, c='blue')
	#plt.show()
		
	return thru_pts
	
#throughput_from_sigmoidfit_coeffs(-40, 0, 10, 20/1, np.linspace(-10,10,100))

## Bayesian linear regression
RegressionResult = namedtuple('RegressionResult',
                                    ['params', 'predict'])

# Assume known fixed prediction error sigma_y
# Prior on w is N(mu_w, sigma_w)
def bayes_lin_reg(X_in, y, mu_w, sigma_w, sigma_y):
    X = np.hstack([np.ones((X_in.shape[0], 1)), X_in]) # add fixed feature for offset; n x d+1

    s = (1.0/sigma_y**2)
    sigma_w_post = np.linalg.inv(np.linalg.inv(sigma_w) + s * X.T @ X)
    mu_w_post = sigma_w_post @ (np.linalg.inv(sigma_w) @ mu_w + s * (X.T @ y))

    def pred(X_in, return_std=False):
        X = np.hstack([np.ones((X_in.shape[0], 1)), X_in])
        if return_std:
            # Instead of einsum we could write :  sigma_y + np.sum(X @ sigma_w_post.* X, axis=1)
            return np.einsum('xd, da -> x', X, mu_w_post), \
                    np.sqrt(sigma_y**2 + np.einsum('xd,de,xe -> x', X, sigma_w_post, X))
        else:
            return X @ mu_w_post
    return RegressionResult((mu_w_post, sigma_w_post), pred)

"""
### linear data
n = 100 ## Play around with n
sigma = 0.2
X, y = data_gen.regression_data_gen(1, 1, lambda x: x/10 + 1, n, sigma, -10, 10)
brr0 = bayes_lin_reg(X, y, np.zeros((2, 1)), np.eye(2), sigma)
print(brr0.params[0], brr0.params[1])
utils.plot_data_and_fit_with_stdev(X, y, brr0.predict, -20, 20)


### Sinusoidal data and Polynomial kernel
pk = 7
n = 10
dim = pk+1
sigma = .2
X_orig, y = data_gen.regression_data_gen(1, 1, lambda x: np.sin(x), n, sigma, -10, 10)
X = utils.polynomial_features(X_orig, pk)
blr2 = bayes_lin_reg(X, y, np.zeros((dim, 1)), 10*np.eye(dim), sigma)
def predict(X_in, return_std=False):
    X = utils.polynomial_features(X_in, pk)
    return blr2.predict(X, return_std)
utils.plot_data_and_fit_with_stdev(X_orig, y, predict, -10, 10)
"""