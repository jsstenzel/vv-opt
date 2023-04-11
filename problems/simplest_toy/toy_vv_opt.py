import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt

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
    
#we'll assume a uniform prior on theta, from 0 to 1
def Ptheta_rvs(s=1):
    #sample from the prior probability of theta
    return scipy.stats.norm.rvs(size=s)
    
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

if False:
    #plot U across the simulated evidence
    d = 0.5
    u, u_list = U_brute_varH(d, eta, H, Gamma, Ptheta_rvs, Ptheta_pdf, N=100, burnin=0, lag=1, verbose=True)
    print(u)
    plt.plot(u_list)
    plt.show()
    uncertainty_prop_plot(u_list)

if False:
    #find optimal d
    u_list = []
    d_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    for d in d_list:
        print("d =",d,'...',flush=True)
        u, _ = U_brute_varH(d, eta, H, Gamma, Ptheta_rvs, Ptheta_pdf, N=300, burnin=0, lag=1)
        u_list.append(u)
        
    plt.plot(d_list, u_list)
    plt.show()
    
if False:
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
if False:
    #plot U across the simulated evidence
    d = 0.5
    u, ulist = U_reloop_varH(d, eta, H, empty, Ptheta_rvs, Ptheta_pdf, n1=1000, n2=500, burnin=0, lag=1) 
    print(u)
    uncertainty_prop_plot(ulist)
    
if True:
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