#pdf estimation methods, including kernel density estimations

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


def multivar_likelihood_kernel(d, exp_fn, prior_rvs, n1):
	thetas = prior_rvs(n1)
	Y1_list = [exp_fn(theta, d) for theta in thetas]
	
	y_theta_values = np.array([np.concatenate([thetas[i],Y1_list[i]]) for i,_ in enumerate(thetas)])
	likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)

	return likelihood_kernel, Y1_list
	
	
#if you want to set the ythetas manually, use this
def general_likelihood_kernel(*params, bw_method='scott'):
	#first ensure all lists have equal length
	if False in [len(i) == len(params[0]) for i in params]:
		print("general_likelihood_kernel failure: all inputs must have same length,", [len(i) for i in params])
		
	#then put lists together
	y_theta_values = np.array([np.concatenate([param[i] for param in params]) for i,_ in enumerate(params[0])])

	likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T, bw_method=bw_method)			
	return likelihood_kernel, y_theta_values

#likelihood is expected to be the kde generated from sample
#this assumes 3 elements in 7 - i can generalize this easily later
def kde_plot(likelihood, sample, plotStyle='separate', center=[], res=1000, ynames=None, plot_xyz=[]):
	data = np.transpose(sample)
	
	if center == []:
		#we have to evaluate at a point in the point in the domain of the pdf
		#so that we can see all the 2d plots at that point
		center = [np.mean(idata) for idata in data]
		
	if ynames == None:
		ynames = ["Y"+str(i) for i,_ in enumerate(data)]
	
	if plotStyle == '3d':
		if plot_xyz==[]:
			plot_xyz=[0,1,2]
		#3d plot
		xdata = [y[plot_xyz[0]] for y in sample]
		ydata = [y[plot_xyz[1]] for y in sample]
		zdata = [y[plot_xyz[2]] for y in sample]
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		ax.scatter3D(xdata, ydata, zdata, c='forestgreen', alpha=0.25)
		plt.show()
	elif plotStyle == 'separate':
		if plot_xyz==[]:
			plot_xyz=range(len(data))
		for i,xdata in enumerate(data):
			if i in plot_xyz:
				xs = np.linspace(min(xdata), max(xdata), res)
				#probs = [likelihood(center but put x at center[i])) for x in xs]
				probs = [likelihood([x if j==i else center[j] for j,_ in enumerate(data)]) for x in xs]
				#plt.plot(xs, [likelihood([x,np.mean(ydata),np.mean(zdata)]) for x in xs])
				plt.axvline(center[i], c='gray')
				plt.plot(xs, probs)
				plt.xlabel(ynames[i])
				plt.show()
	elif plotStyle == 'together':
		if plot_xyz==[]:
			plot_xyz=range(len(data))
		fig, axs = plt.subplots(nrows=1, ncols=len(plot_xyz), sharey='all')
		j=0
		for i,xdata in enumerate(data):
			if i in plot_xyz:
				xs = np.linspace(min(xdata), max(xdata), res)
				#probs = [likelihood(center but put x at center[i])) for x in xs]
				probs = [likelihood([x if j==i else center[j] for j,_ in enumerate(data)]) for x in xs]
				#plt.plot(xs, [likelihood([x,np.mean(ydata),np.mean(zdata)]) for x in xs])
				axs[j].axvline(center[i], c='gray')
				axs[j].plot(xs, probs)
				axs[j].set_xlabel(ynames[i])
				j+=1
		
		plt.subplots_adjust(left=0.1,
							bottom=0.1,
							right=0.9,
							top=0.9,
							wspace=0.0,
							hspace=0.0)
		plt.show()
	else:
		print("kde_plot plotStyle not recognized")
		
		
	
#this function takes in a theta and d,
#runs eta(theta,d) many times,
#grabs the histograms of y of those evaluations
#and coarsely estimates pdf(y) from that
#(This kinda assumes that theres no covariance in the data...)
def eta_histogram_pdf(y, theta, d, eta):
	0
	
#this function assumes that there will be a gaussian form of the likelihood.
#runs eta(theta,d) many times,
#runs np.mean and np.cov on each dimension,
#and calculates returns pdf(y) for the multivariate gaussian
def eta_multigaussian_logpdf(y, theta, d, eta, n_pde=1000):
	ysample = [eta(theta, d) for _ in range(n_pde)]
	means = np.mean(ysample, axis=0)
	
	conditioned_ysample = [[yi/means[i] - 1.0 for i,yi in enumerate(yy)] for yy in ysample]
	conditioned_y = [yi/means[i] - 1.0 for i,yi in enumerate(y)]
	
	covs = np.cov(conditioned_ysample, rowvar=False, ddof=1)
	#for i in len(covs):
	#	covs[i][i] += 1e-6
	return scipy.stats.multivariate_normal.logpdf(conditioned_y, mean=[0 for _ in means], cov=covs)