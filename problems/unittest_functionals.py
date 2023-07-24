import sys
import os
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../..')
from problem_definition import *
from fp_verification.fp_experiment_models import quantum_efficiency_exp


#basics
def dont():
	0


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

qe_mean = Functional([(400,.25),(500,.45),(650,.75),(900,.45),(975,.05)])
qe_mean.set_xlim(350,975)
qe_mean.set_ylim(0,1)
qe_mean.spline_interp(3)
qe_mean.plot(linecolor='k', pointcolor='k', show=False)
plt.show()

#Get average
print(qe_mean.idx(0))
list_wavelength_qe = qe_mean.to_array(10)
print(list_wavelength_qe)
qes = [wave_qe[1] for wave_qe in list_wavelength_qe]
avg_qe = np.mean(qes)
print(avg_qe)


#S_pd = Functional([(200,.12),(280,.1),(300,.125),(400,.185),(500,.25),(600,.32),(633,.33),(800,.45),(900,.4875),(930,.5),(1000,.45),(1100,.15)]) #nm, A/W
S_pd = Functional([(200,.12),(280,.1),(300,.125),(400,.185),(633,.33),(900,.4875),(930,.5),(1000,.45),(1100,.15)]) #nm, A/W
S_pd.set_xlim(350,975)
S_pd.set_ylim(0,1)
S_pd.spline_interp(3)
S_pd.plot()

S_pd_new = noise_to_functional(S_pd, .05)
S_pd_new.plot()



############################
S_pd = Functional([(200,.12),(280,.1),(300,.125),(400,.185),(633,.33),(930,.5),(1000,.45),(1100,.15)]) #nm, A/W
S_pd.set_xlim(350,975)
S_pd.set_ylim(0,1)
S_pd.spline_interp(3)
fp_x_defs = {
				"nx": 2048,
				"ny": 2048,
				"sigma_dc": .5, #e-/s #WAG for now
				"S_pd": S_pd,   #representative photodiode, Hamamatsu S1337-1010BQ
				"S_pd_err": .01  #mA/W
			}
			
S = quantum_efficiency_exp(qe_mean, 1.1, 2.0, 5, 1, .01, fp_x_defs)
S.plot()