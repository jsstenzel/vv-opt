import sys
import os
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


#This sketch is a little optimizer to find a beta distribution with certain properties
#Namely, what properties a,b of the beta distribution give cdf(0.05)=0.95, with a large variance?
#I'm doing this for the FRD prior

#First, define a funciton that is minimized when the two critera are met:
#1. cdf(0.05)=0.95
#2. mean=0.025
mean = 0.025

def min_func(params):
	a = params[0]
	b = (a/mean) - a
	if a<0 or b<0:
		return 1000
	#mean = a/(a+b)
	var = a*b / ( (a+b)*(a+b)*(a+b+1) )

	Prob = scipy.stats.beta.cdf(0.05, a, b, loc=0, scale=1)
	
	criterion1 = (0.95 - Prob)**2 #maximum criterion
	
	return criterion1
	
#solve
res = scipy.optimize.minimize(min_func, [0.5])
aa = res.x[0]
bb = (aa/mean) - aa#res.x[1]
success = res.success
if not success:
	print("Fail")
	sys.exit()

#check
print("a:",aa,"b:",bb)
print("CDF(0.05,a,b)=",scipy.stats.beta.cdf(0.05, aa, bb))
print("Var[Beta(a,b)]=",aa*bb / ( (aa+bb)*(aa+bb)*(aa+bb+1) ))
print("E[Beta(a,b)]=",aa/(aa+bb),flush=True)

xplot = np.linspace(0,1,10000)
cdfplot = [scipy.stats.beta.cdf(x, aa, bb) for x in xplot]
betaplot = [scipy.stats.beta.pdf(x, aa, bb) for x in xplot]
plt.plot(xplot,cdfplot,c='steelblue')
plt.plot(xplot,betaplot,c='green')
plt.xlim(0,1)
plt.title("Plot of the PDF and CDF of the solution")
plt.show()