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

#likelihood is expected to be the kde generated from sample
#this assumes 3 elements in 7 - i can generalize this easily later
def kde_plot(likelihood, sample, plotStyle='separate', center=None):
	data = np.transpose(sample)
	
	if center == None:
		#we have to evaluate at a point in the point in the domain of the pdf
		#so that we can see all the 2d plots at that point
		center = [np.mean(idata) for idata in data]
	
	if plotStyle == '3d':
		#3d plot
		xdata = [y[0] for y in sample]
		ydata = [y[1] for y in sample]
		zdata = [y[2] for y in sample]
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')
	
	if plotStyle == 'separate':
		for i,xdata in enumerate(data):
			xs = np.linspace(min(xdata), max(xdata), len(xdata))
			#probs = [likelihood(center but put x at center[i])) for x in xs]
			probs = [likelihood([x if j==i else center[j] for j,_ in enumerate(data)]) for x in xs]
			#plt.plot(xs, [likelihood([x,np.mean(ydata),np.mean(zdata)]) for x in xs])
			plt.plot(xs, probs)
			plt.show()

	if plotStyle == 'together':
		fig, axs = plt.subplots(nrows=1, ncols=len(data), sharey='all')
		for i,xdata in enumerate(data):
			xs = np.linspace(min(xdata), max(xdata), len(xdata))
			#probs = [likelihood(center but put x at center[i])) for x in xs]
			probs = [likelihood([x if j==i else center[j] for j,_ in enumerate(data)]) for x in xs]
			#plt.plot(xs, [likelihood([x,np.mean(ydata),np.mean(zdata)]) for x in xs])
			axs[i].plot(xs, probs)
			axs[i].set_xlabel("Y"+str(i))
		
		plt.subplots_adjust(left=0.1,
							bottom=0.1,
							right=0.9,
							top=0.9,
							wspace=0.0,
							hspace=0.0)
		plt.show()