import sys
import os
import scipy.stats
import numpy as np

sys.path.append('../..')
from problem_definition import *


#basics
def dont():
	0


#unit testing
theta_simple = [ 
	("theta1", ["funct_splines", [[(350,.25),(1050,.3)],
									1,
									.25,
									.25,
									[0.0,1.0],
									5]])
	]

theta_representative = [ 
	("theta1", ["funct_splines", [[(350,.25),(500,.7),(800,.5),(1050,.3)],
									3,
									.1,
									.1,
									[0.0,1.0],
									1]])
	]


test = ProblemDefinition(dont,dont,dont,theta_simple,[],[],[])

print(test.prior_rvs(1))