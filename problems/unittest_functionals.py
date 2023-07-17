import sys
import os
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../..')
from problem_definition import *


#basics
def dont():
	0

"""
#test Functional class
qe1 = Functional([(350,.25),(1050,.3)])
print(qe1.idx(1))
print(qe1.f(400))

qe2 = Functional([(350,.25),(500,.7),(800,.5),(1050,.3)])
qe2.set_ylim(0,1)
print(qe2.idx(1))
print(qe2.f(400))
qe2.plot()
qe2.spline_interp(3)
qe2.plot()
"""


#unit testing
theta_simple = [ 
	("theta1", ["funct_splines", [[(350,.5),(1050,.5)],
									1,
									.25,
									[350,1050],
									[0.0,1.0],]])
	]
	
test = ProblemDefinition(dont,dont,dont,theta_simple,[],[],[])
funct_sample = test.prior_rvs(1)[0]
funct_sample.plot()

########
theta_representative = [ 
		("unf", ["uniform", [-2,-1]]),
		("red-depletedNIR", ["funct_splines", [[(400,.25),(500,.45),(650,.75),(900,.45),(975,.05)],
									3,
									.05,
									[350,975], #LLAMAS spectral range
									[0.0,1.0],
								 ]])
	]


test = ProblemDefinition(dont,dont,dont,theta_representative,[],[],[])
samples = test.prior_rvs(20)
print(funct_sample)
for sample in samples:
	funct = sample[1]
	funct.plot(show=False)
	
#plot mean
qe_mean = Functional([(400,.25),(500,.45),(650,.75),(900,.45),(975,.05)])
qe_mean.set_xlim(350,975)
qe_mean.set_ylim(0,1)
qe_mean.spline_interp(3)
qe_mean.plot(linecolor='r', pointcolor='r', show=False)
	
plt.show()