import sys
import os
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('..')

#This sketch is intended to estimate the convergence of sample probability
#sample probability is what im referring to as an estimator of a probability:
#P(x < a) for x~X =~ (1/N)SUM(xi < a) for xi~X

#just using gaussians for this for now

def sample_prob(a, N): 
	true_prob = scipy.stats.norm.cdf(a, loc=0, scale=1)
	
	samples = scipy.stats.norm.rvs(size=N, loc=0, scale=1)
	num_true = 0
	for xi in samples:
		num_true += int(xi <= a)
	est_prob = num_true / N
	
	return true_prob, est_prob
	

a = scipy.stats.norm.rvs(size=1, loc=0, scale=1)[0]
print(sample_prob(a, 10))
print(sample_prob(a, 10**2))
print(sample_prob(a, 10**3))
print(sample_prob(a, 10**4))
print(sample_prob(a, 10**5))
print(sample_prob(a, 10**6))

ns = [1e0,1e1,1e2,1e3,1e4,1e5,1e6]
ns += [3e0,3e1,3e2,3e3,3e4,3e5,3e6]
ns.sort()

estimates=[sample_prob(a, int(np.floor(n))) for n in ns]
diffs = [abs(est[0] - est[1]) for est in estimates]

plt.xscale('log')
plt.plot(ns, estimates)
plt.show()

plt.xscale('log')
plt.plot(ns, diffs)
plt.show()