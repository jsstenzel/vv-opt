import sys
import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

sys.path.append('..')
from totalOrderMultiIndexSet import *

"""
Taking a similar approach as detailed in Huan & Marzouk 2013
But not using dimension-adaptive sparse quadrature
"""
def likelihood_pce(problem, p):
	stochastic_dim = problem.dim_theta + problem.dim_d
	ny = problem.dim_y
	total_multindex = totalOrderMultiIndices(stochastic_dim, p)
	
	#construct the psis for each multindex ii
	def Psi(x, ii):
		psi = 1
		for i in ii: #each multindex defines the index of each parameter in stochastic_dim
			if i < problem.dim_theta :
				#assume it's Gaussian, over the whole domain, use Hermite
				psi*=Hermite(x,i)
			else:
				#assume its uniform, from 0 to 1, use Legendre
				#(transform the domain of theta elsewhere)
				psi*=Legendre(x,i)
		
		return psi
	
	total_G = []
	for ii in total_multindex:
		G_ii = []
		psi_sq_p
		
		for c in range(ny):
			#Solve the integral in Equation 21
			
			G_ii_c = integral / psi_sq_p

	
	return total_G
	
def sample_likelihood_pce(problem, total_G):
	stochastic_dim = problem.dim_theta + problem.dim_d
	total_multindex = totalOrderMultiIndices(stochastic_dim, p)

def eval_likelihood_pce(problem, total_G):
	stochastic_dim = problem.dim_theta + problem.dim_d
	total_multindex = totalOrderMultiIndices(stochastic_dim, p)

def Legendre(x, a):
	c = np.zeros(a+1)
	c[a] = 1
	psi_a = np.polynomial.legendre.legval(x, c)
	return psi_a

def Hermite(x, a):
	c = np.zeros(a+1)
	c[a] = 1
	psi_a = np.polynomial.hermite_e.hermeval(x, c)
	return psi_a


################################################################









def K_lognormal_fn(Z, mu=1, b=[.2,.4,.2,.4,.2]):
	b_mat = np.array(b).T
	Z_mat = np.array(Z)
	power = float(b_mat @ Z_mat + mu)
	return np.exp(power)
	
def sample_K_lognormal(mu=1, b=[.2,.4,.2,.4,.2]):
	Z = scipy.stats.norm.rvs(size=len(b))
	return K_lognormal_fn(Z, mu, b)
	
def plot_hist_K_samples(n, showplot=True):
	samples = [sample_K_lognormal() for i in range(n)]
	uncertainty_prop_plot(samples,xlab="K", c='maroon', title="True samples of K",showFig=showplot)

#plot_hist_K_samples(10000)
	
#Use Hermite polynomials to approximate the function of random variables
#each of the n random variables Z have a truncation in the multi-index alpha
#using the recommended truncation scheme, with the provided multi index enumeration function
def poly_approx_K(p, mu=1, b=[.2,.4,.2,.4,.2]):
	n = len(b)
	#First, enumerate the multi indices
	total_multindex = totalOrderMultiIndices(n, p)

	#Now we have a list of all of the different alphas we want to include
	#Next, iterate over the alphas in total_multi_index
	#For each alpha, calculate the polynomial coefficient, and add to a list
	total_c_alpha = []
	for alpha in total_multindex:
		c_alpha = np.exp(mu)
		for i,ai in enumerate(alpha):
			c_alpha *= (b[i]**ai / math.factorial(ai)) * np.exp(b[i]**2 / 2)
		
		total_c_alpha.append(c_alpha)
	
	#Lastly, return those coefficients, so you can use them to evaluate the polynomial approx elsewhere
	return total_c_alpha, total_multindex
	
#print(poly_approx_K(1))
#print(poly_approx_K(2))
	
#Coded my own function because I didn't like the provided one
#Just calculates He_a(x)

	
def eval_poly_approx_K(total_c_alpha, total_multindex, Z=[]):
	n = len(total_multindex[0])
	p = total_multindex[-1][0] #gross but idc it works
	if len(Z) == 0:
		Z = scipy.stats.norm.rvs(size=n)

	#Then, with all of the c_alpha in hand, we want to evaluate the polynomial approximate
	evaluation = 0
	for j,c_alpha in enumerate(total_c_alpha):
		alpha = total_multindex[j]
		#print(alpha)
		multi_Hermite = np.prod([Hermite(Z[i], ai) for i,ai in enumerate(alpha)])
		evaluation += c_alpha * multi_Hermite
		
	return evaluation
	
#c,mi  = poly_approx_K(1)
#print(eval_poly_approx_K(c, mi, Z=[0.1,0.1,0.1,0.1,0.1]),flush=True)

	
def plot_hist_K_approx(n, p, showplot=True):
	#First, get the coefficients
	c,mi  = poly_approx_K(p)
	
	#Then sample
	samples = [eval_poly_approx_K(c,mi) for i in range(n)]
	uncertainty_prop_plot(samples,xlab="K", title="Approximate samples of K (p="+str(p)+")",showFig=showplot)