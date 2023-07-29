import sys
import os
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import math

sys.path.append('../..')
from problem_definition import *
from fp_verification.fp_experiment_models import quantum_efficiency_exp


#basics
def dont():
	0


#test Functional class
if False:
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

if False:
	test = ProblemDefinition(dont,dont,dont,theta_simple,[],[],[])
	funct_sample = test.prior_rvs(1)[0]
	funct_sample.plot()

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
if False:
	print(qe_mean.idx(0))
	list_wavelength_qe = qe_mean.to_array(10)
	print(list_wavelength_qe)
	qes = [wave_qe[1] for wave_qe in list_wavelength_qe]
	avg_qe = np.mean(qes)
	print(avg_qe)

if False:
	#S_pd = Functional([(200,.12),(280,.1),(300,.125),(400,.185),(500,.25),(600,.32),(633,.33),(800,.45),(900,.4875),(930,.5),(1000,.45),(1100,.15)]) #nm, A/W
	S_pd = Functional([(200,.12),(280,.1),(300,.125),(400,.185),(633,.33),(900,.4875),(930,.5),(1000,.45),(1100,.15)]) #nm, A/W
	S_pd.set_xlim(350,975)
	S_pd.set_ylim(0,1)
	S_pd.spline_interp(3)
	S_pd.plot()

	S_pd_new = noise_to_functional(S_pd, .05)
	S_pd_new.plot()



############################
if False:
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
	
	
#unit test for multivariate likelihood kernel
if True:
	"""
	def calc_likelihood_kernel(d, exp_fn, p_theta_rvs, n1=1000, c='r', showplot=True):
		thetas = p_theta_rvs(n1)
		Y1_list = [exp_fn(theta, d) for theta in thetas]
		
		y_theta_values = np.vstack([Y1_list, thetas])
		likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values)

		#lets plot this real quick
		if showplot:
			sns.kdeplot(Y1_list, thetas, color=c, shade=True, cmap="Reds", shade_lowest=False)
			plt.show()

		return likelihood_kernel, Y1_list
	"""
	from problems.fp_verification.fp_vv_opt import *
	test_d = [
				20,   #t_gain
				30,   #I_gain
				5,    #n_meas_rn
				6,    #d_num
				7200, #d_max
				1     #d_pow
			 ]

	thetas = fp.prior_rvs(100000)
	#print("thetas")
	#print(*thetas,sep='\n')
	
	ys = [fp.eta(theta, test_d) for theta in thetas]
	#print("\nys")
	#print(*ys,sep='\n')
	
	y_theta_values = np.array([np.concatenate([thetas[i],ys[i]]) for i,_ in enumerate(thetas)])
	#print("\ny_theta_values")
	#print(y_theta_values)
		
	likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)
	#needs transpose apparently for multivariate dataset. Tricky...
	"""
	#lets test this real quick: passing in same thetas and new ys should result in high probs
	new_ys = [fp.eta(theta, test_d) for theta in thetas] #different due to exptl. noise
	new_ythetas = np.array([np.concatenate([thetas[i],new_ys[i]]) for i,_ in enumerate(thetas)])
	result = likelihood_kernel(new_ythetas.T)
	print("\nhigh likelihood results")
	print(result)
	
	#resampling and checking likelihood should produce okay results?
	_thetas = fp.prior_rvs(10)
	_ys = [fp.eta(theta, test_d) for theta in _thetas] #different due to exptl. noise
	_ythetas = np.array([np.concatenate([_thetas[i],_ys[i]]) for i,_ in enumerate(thetas)])
	result = likelihood_kernel(_ythetas.T)
	print("\nok likelihood results")
	print(result)
	
	#now lets test it with a big bias on a theta?
	bad_ythetas = np.array([ytheta+[.1,0,0,0,0,0] for ytheta in new_ythetas])
	result = likelihood_kernel(bad_ythetas.T)
	print("\nlow likelihood results")
	print(result)
	"""
	#check a single
	theta = fp.prior_rvs(1)
	y = fp.eta(theta, test_d)
	ytheta = np.concatenate([theta, y])
	print("\n ok likelihood result")
	print(ytheta)
	print(likelihood_kernel(ytheta))
	print(likelihood_kernel(ytheta.T))
	theta_pdfs = fp.prior_pdf_unnorm(theta)
	print("theta pdfs")
	print(theta_pdfs)
	print(np.prod(theta_pdfs))
	
	#test with theta and y_nominal
	theta_nominal = [1.1, 2.5, .001]
	theta = fp.prior_rvs(1)
	y_nominal = fp.eta(theta_nominal, test_d)
	ytheta = np.concatenate([theta, y_nominal])
	print("\nlikelihood result of unrelated theta and y")
	print(ytheta)
	print(likelihood_kernel(ytheta))

#gain experiment testing
if False:
	from scipy.stats import norm, poisson

	#Assumptions:
	#temperature is certain
	#w is certain
	#experiment runtimes are certain
	_temp= -90+273.15 #K
	_k = 1.380649e-23 #J / K
	_c = 299792458 #m / s
	_e0 = 22100 #eV
	_m0 = 108.9049856 * 1.66054e-27 #cd109 atomic mass in kg
	_An = 6.02214076e23 #avogadros constant
	xlist = [
			#general
			("nx", 2048),
			("ny", 2048),
			("sigma_dc", .5), #e-/s #Estimate based on a consistency test performed on SN20006
			("mu_stray", .1), #e-/s #WAG for now
			("sigma_stray", .005), #WAG for now
			#gain
			("P_signal", 0.90), #Prob. of correctly identifying signal as event #WAG for now
			("P_noise", 0.01), #Prob. of incorrectly identifying noise/interference as event #WAG for now
			("T_ccd", _temp), #K
			("E0", _e0), #22.1 keV Cd-109 emission line
			("sigma_E", math.sqrt((_m0 * _c**2) / (_k*_temp*_e0**2))), #1/eV^2
			("w", 3.66 + 0.000615*(300-_temp)), #eV/e- #this has uncertainty. nuisance parameter?
			("activity_cd109", 5e-6), #Ci #radioactivity of sample
			("grade_size", 3), #3x3 event grade sizes
			("t_gain_setup", 1200), #WAG
			("t_gain_buffer", 5), #WAG
			#rn
			("t_rn", .1), #100 ms exposure
			("t_rn_buffer", 5), #WAG
			#dc
			("t_0", 0.1), #100ms baseline exposure assumed
			("t_dc_buffer", 5), #WAG
			#qoi
			("tau", 1800),
			#cost
			("testbed_setup", 1800), #WAG
			("C_engineer", 0.00694444444) #WAG $/s, from $25/hr
		 ]
	_x = dict(xlist)
	#define parameters:
	P_signal = _x["P_signal"] #probability of a signal event being correctly identified as an event; 1-P_signal is false negative rate
	P_noise = _x["P_noise"] #probability of noise being incorrectly identified as an event; P_noise is false positive rate
	sigma_dc = _x["sigma_dc"]
	T_ccd = _x["T_ccd"]
	E0 = _x["E0"]
	sigma_E = _x["sigma_E"]
	w = _x["w"]
	activity_cd109 = _x["activity_cd109"] #look at cd-109 sample
	nx = _x["nx"]
	ny = _x["ny"]
	grade_size = _x["grade_size"] #we are considering 3x3 event grade sizes, based on the physics of the experiment
	
	#define design variables
	t=20
	I=30
	#other args
	gain = 1.1
	rn = 2.5
	dc = 0.001
	
	#Sample from Poisson distribution of # particle events
	#no, doing this non-randomly at first
	#number of isotopes after time t: 
	#N_t = N0_cd109 * math.exp(-rate_cd109 * t)
	#number of decays after time t: activity = rate_cd109 * N0_cd109
	Activity_s = activity_cd109 * 3.7e10 #convert to decay/s
	N_decays = Activity_s * t
	print("Activity_s", Activity_s)
	print("N_decays", N_decays)
	
	#calculate derived values
	grade_area = grade_size**2
	n_signal = N_decays * P_signal
	n_noise = (nx*ny/grade_area) * P_noise
	n = n_signal + n_noise
	p = n_signal / n
	print("n_signal",n_signal)
	print("n_noise",n_noise)
	print("p", p)
	
	#calculate mixture distribution parameters, X = pS + (1-p)N
	mu_x = p*E0*gain/w
	var_x = p * (sigma_E/w)**2 + rn**2 + sigma_dc**2 + p*(1-p)*(E0/w)**2
	#Im missing I here
	sigma_x = math.sqrt(var_x)
	print("sigma_x",sigma_x)
	stddev = sigma_x/math.sqrt(I)
	print("stddev",stddev)
	
	random = norm.rvs(scale = stddev)
	y = mu_x + random
	print("mu_x",mu_x)
	print("random",random)