import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv

sys.path.append('../..')
from obed.obed import *
from uq.uncertainty_propagation import uncertainty_prop_plot
from obed.convergence_solver import *

def eta(_theta, _d):
    random = scipy.stats.norm.rvs(loc=0, scale=_d, size=1)[0]
    return _theta + random
    
def H(_theta):
    return 10*_theta**2
    
def Gamma(_d):
    return 10/_d
    
def Ptheta_rvs(s=1):
    #sample from the prior probability of theta
    return scipy.stats.norm.rvs(size=s, loc=0.5, scale=0.5)

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
    d_list = [.05,.1,.15,.2,.25,.3,.35,.4,.45,.5]
    G_list = []
    for d in d_list:
        print("d =",d,'...',flush=True)
        u, _ = U_reloop_varH(d, eta, H, empty, Ptheta_rvs, Ptheta_pdf, n1=3000, n2=1000, burnin=0, lag=1)
        u_list.append(u)
        G_list.append(Gamma(d))
        
    with open('reloop_u_plot_pareto.csv', 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(u_list)
        
    plt.plot(d_list, u_list, 'g')
    plt.plot(d_list, G_list, 'r')
    plt.xlabel("d")
    plt.ylabel("u, Γ")
    plt.show()
    plt.plot(G_list,u_list)
    plt.xlabel("Γ")
    plt.ylabel("u")
    plt.show()
    
    

#Calls for some of the convergence functions:
#trace_scipy_kde(d=0.5, eta, H, Gamma, Ptheta_rvs, Ptheta_pdf, nmax=10000)

#trace_loop2(d=0.5, y=0.6, exp_fn, H, Gamma, Ptheta_rvs, Ptheta_pdf, burnin=0, lag=1, min=100, doLog=True, goPlot=True)

#analyze_trace_loop2("loop2.csv")

#loop3_convergence(d=0.5, exp_fn, H, Gamma, Ptheta_rvs, Ptheta_pdf, calc_likelihood_kernel, doLog=True, doPlot=True)
