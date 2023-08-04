#pdf estimation methods, including kernel density estimations

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns


		

def multivar_likelihood_kernel(d, exp_fn, prior_rvs, n1):
	thetas = prior_rvs(n1)
	Y1_list = [exp_fn(theta, d) for theta in thetas]
	
	y_theta_values = np.array([np.concatenate([thetas[i],Y1_list[i]]) for i,_ in enumerate(thetas)])
	likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)

	return likelihood_kernel, Y1_list
	
	
#if you want to set the ythetas manually, use this
def general_likelihood_kernel(*params):
	#first ensure all lists have equal length
	if False in [len(i) == len(params[0]) for i in params]:
		print("general_likelihood_kernel failure: all inputs must have same length,", [len(i) for i in params])
		
	#then put lists together
	y_theta_values = np.array([np.concatenate([param[i] for param in params]) for i,_ in enumerate(params[0])])

	likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)			
	return likelihood_kernel

#likelihood is the kde generated from sample
#this assumes 3 elements in 7 - i can generalize this easily later
def kde_plot(likelihood, sample):
	xdata = [y[0] for y in sample]
	ydata = [y[1] for y in sample]
	zdata = [y[2] for y in sample]
	
	#3d plot
	#fig = plt.figure()
	#ax = plt.axes(projection='3d')
	#ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')
	
	xs = np.linspace(min(xdata), max(xdata), len(xdata))
	plt.plot(xs, [likelihood([x,np.mean(ydata),np.mean(zdata)]) for x in xs])
	plt.show()
	
	ys = np.linspace(min(ydata), max(ydata), len(ydata))
	plt.plot(ys, [likelihood([np.mean(xdata),y,np.mean(zdata)]) for y in ys])
	plt.show()
	
	zs = np.linspace(min(zdata), max(zdata), len(zdata))
	plt.plot(zs, [likelihood([np.mean(xdata),np.mean(ydata),z]) for z in zs])
	plt.show()