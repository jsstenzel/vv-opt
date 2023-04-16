import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sys.path.append('../..')
from obed.obed import *
from uq.uncertainty_propagation import uncertainty_prop_plot
from obed.convergence_solver import *

def eta(_theta, _d):
    random = scipy.stats.norm.rvs(loc=0, scale=_d, size=1)[0]
    return _theta + random
    
def H(_theta):
    return (_theta+1)**2
    
def Gamma(_d):
    return 10/(_d+.01)
    
def Ptheta_rvs(s=1):
    #sample from the prior probability of theta
    return scipy.stats.norm.rvs(size=s, loc=0.5, scale=0.5)
    
#def Ptheta_pdf(_theta):
#    #return the prior probability of theta
#    if _theta<0 or _theta>1:
#        return 0
#    else:
#        return 1

def Ptheta_pdf(_theta):
    return scipy.stats.norm.pdf(_theta, loc=0.5, scale=0.5)
    
def empty(whatever):
    return 0

#############################################

def __brute_u_example():
    #plot the trace of H(posterior)
    d = 0.5
    u, H_trace = U_brute_varH(d, eta, H, Gamma, Ptheta_rvs, Ptheta_pdf, N=1000, burnin=0, lag=1)
    print(u)
    plt.plot(H_trace)
    plt.show()
    uncertainty_prop_plot(H_trace)

def __brute_u_find_d():
    #find optimal d
    u_list = []
    d_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    for d in d_list:
        print("d =",d,'...',flush=True)
        u, _ = U_brute_varH(d, eta, H, Gamma, Ptheta_rvs, Ptheta_pdf, N=300, burnin=0, lag=1)
        u_list.append(u)
        
    plt.plot(d_list, u_list)
    plt.show()
    

def __brute_u_plot_pareto():
    #find optimal d, plotting out u and Gamma
    u_list = []
    d_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    G_list = []
    for d in d_list:
        print("d =",d,'...',flush=True)
        u, _ = U_brute_varH(d, eta, H, empty, Ptheta_rvs, Ptheta_pdf, N=1000, burnin=0, lag=1, verbose=True)
        u_list.append(u)
        G_list.append(Gamma(d))
        
    plt.plot(d_list, u_list, 'g')
    plt.plot(d_list, G_list, 'r')
    plt.xlabel("d")
    plt.ylabel("u, Γ")
    plt.show()
    plt.plot(G_list,u_list)
    plt.xlabel("Γ")
    plt.ylabel("u")
    plt.show()
    
###################################################################

def __reloop_u_example():
    #plot the trace of H(posterior)
    d = 0.5
    u, ulist = U_reloop_varH(d, eta, H, empty, Ptheta_rvs, Ptheta_pdf, n1=1000, n2=500, burnin=0, lag=1) 
    print(u)
    uncertainty_prop_plot(ulist)
    

def __reloop_u_plot_pareto():
    #find optimal d, plotting out u and Gamma
    u_list = []
    d_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    G_list = []
    for d in d_list:
        print("d =",d,'...',flush=True)
        u, _ = U_reloop_varH(d, eta, H, empty, Ptheta_rvs, Ptheta_pdf, n1=1000, n2=500, burnin=0, lag=1)
        u_list.append(u)
        G_list.append(Gamma(d))
        
    plt.plot(d_list, u_list, 'g')
    plt.plot(d_list, G_list, 'r')
    plt.xlabel("d")
    plt.ylabel("u, Γ")
    plt.show()
    plt.plot(G_list,u_list)
    plt.xlabel("Γ")
    plt.ylabel("u")
    plt.show()
    
    
#To try to answer how large n1 needs to be, im going to try plotting trace of gaussian kde
def __trace_scipy_kde():
    nmax = 10000
    d = 0.5
    thetas = Ptheta_rvs(100)
    ys = [eta(theta, d) for theta in thetas]
    
    test_thetas = Ptheta_rvs(100)
    test_y_theta = [[eta(theta, d),theta] for theta in test_thetas]
    trace = []
    for n in range(nmax):
        thetas = np.append(thetas, Ptheta_rvs(1)[0])
        ys = np.append(ys, eta(thetas[-1], d))
        
        y_theta_values = np.vstack([ys, thetas])
        likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values)
        
        #grab metrics from this kde
        #i want to trace some kind of change from one step to the next
        test_calcs = [likelihood_kernel.pdf(pair)[0] for pair in test_y_theta]
        
        trace_metric = sum(test_calcs)
        trace.append(trace_metric)
    #print(trace, flush=True)
    plt.plot(trace)
    plt.show()


def __loop3_convergence():
	#Here, I want to figure out how many likelihood samples are needed to make a good model of likelihood fn
	#seeking convergence of T=SUM(likelihood(y_test,theta_test)) across likelihood samples
    
	#Come up with a smarter way of determining the fitness of a kde
	0
	
def __loop2_convergence():
	#How many MCMC steps are needed to generate a good sample of the posterior?
	#basically i want to go until Ui is converged
	p_theta_rvs = Ptheta_rvs
	p_theta_pdf = Ptheta_pdf
	exp_fn = eta
	H_fn = H
	d = 0.5
	y = 0.6#exp_fn(p_theta_rvs(1)[0], d) #whatevs just pick one from likelihood
	burnin=0
	lag=1
	min=100
	
	#assume we have a good likelihood kernel to use in loop3. set it up to reuse here
	n1 = 10000
	likelihood_kernel, _ = calc_likelihood_kernel(d, exp_fn, p_theta_rvs, n1, showplot=False)
	
	#setup loop
	theta_current = p_theta_rvs(1)[0]
	mcmc_trace = []
	U_trace = []
	i = 0
	
	while not is_algorithm_converged(U_trace, min, ci=0.975, closeness=0.99):
		#one MCMC loop step
		#Propose a new value of theta from the prior
		theta_proposal = p_theta_rvs(1)[0]
		#Compute likelihood of the "data" for both thetas - reusing the loop-1 likelihood evaluations
		likelihood_current = likelihood_kernel.evaluate([y,theta_current]) #using kde to estimate pdf
		likelihood_proposal = likelihood_kernel.evaluate([y,theta_proposal]) #calc probability at the data y
		#Compute acceptance ratio
		prior_current = p_theta_pdf(theta_current)
		prior_proposal = p_theta_pdf(theta_proposal)
		p_current = likelihood_current * prior_current
		p_proposal = likelihood_proposal * prior_proposal
		R = p_proposal / p_current
		#Accept our new value of theta according to the acceptance probability R:
		if np.random.random_sample() < R:
			theta_current = theta_proposal
		#Include theta_current in my trace according to rules
		i += 1
		if not( i > burnin and i%lag == 0):
			continue #i.e. run the mcmc part again, dont add corresponding U to the trace
		mcmc_trace.append(theta_current)
		
		#We also want to use min here to make sure that we're not taking a meaningless variance for our utility, with too little mcmc data
		if len(mcmc_trace) <= min:
			continue
			
		#Now, use my posterior samples to calculate H(theta|y,d) samples
		H_theta_posterior = [H_fn(theta) for theta in mcmc_trace]
		
		#Now find the variance of that. Return the inverse as the utility; consider cost separately
		var_H = np.var(H_theta_posterior, ddof=1) #its a sample variance
		U = (1.0/var_H)
		U_trace.append(U)
		#print(len(U_trace), end='\r', flush=True)
		
		#more informative print statement
		"""
		sample_mean = np.mean(U_trace)
		sample_var = np.var(U_trace, ddof=1)
		N_latest = len(U_trace)
		print(N_latest, U, sample_mean, sample_var, scipy.stats.norm.ppf(0.95), scipy.stats.norm.ppf(0.95)*np.sqrt(sample_var/N_latest)/sample_mean, flush=True)
		"""
		
	print("N for convergence is", len(U_trace), ", mean U is", np.mean(U_trace))
	
__loop2_convergence()