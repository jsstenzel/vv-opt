import sys
import numpy as np
import numpy as np
import scipy
import matplotlib.pyplot as plt


def kl_divergence_2gaussians(mu1, sig1, mu2, sig2):
	sratio = sig1**2 / sig2**2
	mratio = (mu1 - mu2)**2 / sig2**2
	return 0.5 * (sratio + mratio - 1 - np.log(sratio))
	
def kl_divergence_2gammas(a1, b1, a2, b2):
	return a2*np.log(b1/b2) - scipy.special.gammaln(a1) + scipy.special.gammaln(a2) + (a1-a2)*scipy.special.digamma(a1) - (b1 - b2)*(a1/b1)

def kl_divergence_2bernoullis(q1, q2):
	return q2 * np.log(q2/q1) + (1 - q1) * ln((1-q2)/(1-q1))
	
def make_ecdf(data):
	x = np.sort(data)
	n = len(x)
	x_u, x_rle = np.unique(x, return_counts=True)
	#y = (np.cumsum(x_rle) - 0.5) / n #a peculiarity of the method?
	y = np.cumsum(x_rle) / n
	
	#y_ecdf = np.arange(1, n + 1) / n
	
	ecdf = scipy.interpolate.interp1d(x_u, y, kind='linear', bounds_error=False, fill_value=(0, 1))

	return ecdf
	
#Perez-Cruz 2008
#TODO rewrite this in a way that is robust to small e
def kl_divergence_1dsamplingestimator(s1, s2):
	np.seterr(all='ignore') 
	samples1 = np.float64(s1)
	samples2 = np.float64(s2)
	
	#get e and n
	dx = np.diff(np.sort(np.unique(samples1)))
	dy = np.diff(np.sort(np.unique(samples2)))
	ex, ey = np.min(dx), np.min(dy)
	e = np.float64(min(ex, ey) * 0.5)
	n = len(samples1) 
	print("e:",e)
	
	#get the empirical cdfs
	P = make_ecdf(samples1)
	Q = make_ecdf(samples2)
	
	#summand = [np.log( (P(x)-P(x-e))/(Q(x)-Q(x-e)) ) for x,y in zip(samples1,samples2)]
	summand = [0]*len(samples1)
	drops = 0
	for i,x in enumerate(np.sort(np.unique(samples1))):
		num = P(x)-P(x-e)
		denom = Q(x)-Q(x-e)
		elem_i = np.log(num/denom)
		if denom <= 0 or elem_i == -np.inf or elem_i == np.inf:
			#print("Found a 0 in denominator at",x)
			drops += 1
			n -= 1 #not sure if this fixes or introduces a bias
			continue
		summand[i] = elem_i
	
	D_est = sum(summand) / n   
	print("number of drops:",drops)
	
	return D_est - 1

"""
def kl_divergence_1dsamplingestimator_pretrained(s1, Q, ey):
	samples1 = np.float64(s1)
	
	#get e and n
	dx = np.diff(np.sort(np.unique(samples1)))
	ex = np.min(dx)
	e = np.float64(min(ex, ey) * 0.5)
	n = len(samples1) 
	
	#get the empirical cdfs
	P = make_ecdf(samples1)
	
	#summand = [np.log( (P(x)-P(x-e))/(Q(x)-Q(x-e)) ) for x,y in zip(samples1,samples2)]
	summand = [0]*len(samples1)
	for i,x in enumerate(samples1):
		num = P(x)-P(x-e)
		denom = Q(x)-Q(x-e)
		if denom <= 0:
			#print("Found a 0 in denominator at",x)
			#n -= 1
			continue
		else:
			elem_i = np.log(num/denom)
			summand[i] = elem_i
	
	D_est = sum(summand) / n   
	
	return D_est - 1
"""

if __name__ == '__main__':  
	#test it. Get true KL divergence of two gaussians:
	mu1 = 1
	sig1 = 4
	mu2 = 2
	sig2 = 3
	true_D = kl_divergence_2gaussians(mu1, sig1, mu2, sig2)
	
	#sample from those gaussians
	num_vals = 10**8
	samples1 = scipy.stats.norm.rvs(size=num_vals, loc=mu1, scale=sig1)
	samples2 = scipy.stats.norm.rvs(size=num_vals, loc=mu2, scale=sig2)
	
	#check the ecdfs
	"""
	P = make_ecdf(samples1)
	Q = make_ecdf(samples2)
	minrange = min(min(samples1),min(samples2))
	maxrange = max(max(samples1),max(samples2))
	cushion = (maxrange-minrange) / 20
	ecdf_range = np.linspace(minrange-cushion,maxrange+cushion,num_vals*3)
	Py = [P(x) for x in ecdf_range]
	Qy = [Q(x) for x in ecdf_range]
	plt.plot(ecdf_range, Py, 'r')
	plt.plot(ecdf_range, Qy, 'b')
	plt.show()
	"""
	
	#estimate
	print("True D_KL:",true_D)
	for n in [10,30,100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000,30000000,100000000]:
		s1 = samples1[:n]
		s2 = samples2[:n]
		D_est = kl_divergence_1dsamplingestimator(s1, s2)
		print("n =",n,"Estimated D_KL:",D_est)
	
	#plot it