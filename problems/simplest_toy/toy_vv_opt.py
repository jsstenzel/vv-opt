import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sys.path.append('../..')
from obed.obed import U_brute_varH, U_reloop_varH
from uq.uncertainty_propagation import uncertainty_prop_plot

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
        u, _ = U_brute_varH(d, eta, H, empty, Ptheta_rvs, Ptheta_pdf, N=1000, burnin=0, lag=1)
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
        
__trace_scipy_kde()


def __loop3_convergence():
	#Here, I want to figure out how many likelihood samples are needed to make a good model of likelihood fn
	#seeking convergence of T=SUM(likelihood(y_test,theta_test)) across likelihood samples
    d = 0.5
	test_pop = 100
	burn_in = 4000 #these looked noisy, we dont care about those early matrics
	
    thetas = Ptheta_rvs(burn_in)
    ys = [eta(theta, d) for theta in thetas]

    test_thetas = Ptheta_rvs(test_pop)
    test_y_theta = [[eta(theta, d),theta] for theta in test_thetas]
    trace = []
    while converged == False:
        thetas = np.append(thetas, Ptheta_rvs(1)[0])
        ys = np.append(ys, eta(thetas[-1], d))
        
        y_theta_values = np.vstack([ys, thetas])
        likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values)
        
        #grab metrics from this kde
        test_calcs = [likelihood_kernel.pdf(pair)[0] for pair in test_y_theta]
        trace_metric = sum(test_calcs)
		
        trace.append(trace_metric)
				
		#Test convergence
		converged = is_algorithm_converged(trace, ci=0.95, closeness=0.95)
	
	
	
	
	
	
	