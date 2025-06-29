import sys
import math
import numpy as np
from scipy.stats import norm, poisson

sys.path.append('../..')
from approx.gaussian_process import *

"""
Full matrix experiment model
"""
def fp_likelihood_fn(theta, d, x, dummy, err=True):
	#define interest params:
	if type(theta) is dict:
		gain = theta["gain"]
		rn = theta["rn"]
		dc = theta["dc"]
	else:
		gain = theta[0]
		rn = theta[1]
		dc = theta[2]
	#define design variables, enforcing discreteness:
	t_gain = d["t_gain"]
	I_gain = d["I_gain"]
	n_meas_rn = d["n_meas_rn"]
	d_num = d["d_num"]
	d_max = d["d_max"]
	d_pow = d["d_pow"]
	#just pass along entire x
	
	y1 = gain_exp(gain, rn, dc, t_gain, I_gain, x, err)
	y2 = read_noise_exp(gain, rn, n_meas_rn, x, err)
	y3 = dark_current_exp(gain, rn, dc, d_num, d_max, d_pow, x, err)
	
	return [y1, y2, y3]
	
"""
Full matrix experiment model w/ quantum efficiency
"""
def fp_qe_likelihood_fn(theta, d, x, err=True):
	#define interest params:
	if type(theta) is dict:
		gain = theta["gain"]
		rn = theta["rn"]
		dc = theta["dc"]
		qe = theta["qe"]
	else:
		gain = theta[0]
		rn = theta[1]
		dc = theta[2]
		qe = theta[3]
	#define design variables, enforcing discreteness:
	t_gain = d["t_gain"]
	I_gain = d["I_gain"]
	n_meas_rn = d["n_meas_rn"]
	d_num = d["d_num"]
	d_max = d["d_max"]
	d_pow = d["d_pow"]
	n_qe = d["n_qe"]
	t_qe = d["t_qe"]
	I_qe = d["I_qe"]
	#just pass along entire x
	
	y1 = gain_exp(gain, rn, dc, t_gain, I_gain, x, err)
	y2 = read_noise_exp(gain, rn, n_mean_rn, x, err)
	y3 = dark_current_exp(gain, rn, dc, d_num, d_max, d_pow, x, err)
	y4 = quantum_efficiency_exp(qe, gain, rn, n_qe, t_qe, I_qe, x, err)
	
	return [y1, y2, y3, y4]


"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
Assumptions:
* All decays of the sample make it to the CCD
* Signal and noise events will never overlap
* w has no uncertainty
"""
def gain_exp(gain, rn, dc, t, I, _x, err=True):
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
	#t length of sample exposure
	#I number of exposures
	
	#Sample from Poisson distribution of # particle events
	#no, doing this non-randomly at first
	#number of isotopes after time t: 
	#N_t = N0_cd109 * math.exp(-rate_cd109 * t)
	#number of decays after time t: activity = rate_cd109 * N0_cd109
	Activity_s = activity_cd109 * 3.7e10 #convert to decay/s
	N_decays = Activity_s * t
	
	#calculate derived values
	grade_area = grade_size**2
	n_signal = N_decays * P_signal
	n_noise = (nx*ny/grade_area) * P_noise
	n = n_signal + n_noise
	p = n_signal / n
	
	#calculate mixture distribution parameters, X = pS + (1-p)N
	mu_x = p*E0*gain/w
	var_x = p * (sigma_E/w)**2 + rn**2 + sigma_dc**2 + p*(1-p)*(E0/w)**2
	sigma_x = math.sqrt(var_x)
	stddev = sigma_x/math.sqrt(I)
	
	random = norm.rvs(scale = stddev)
	if err:
		y = mu_x + random
	else:
		y = mu_x
	return y
	

"""
theta: [0] gain [1] read noise [2] dark current
d: [0] number exposures
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def read_noise_exp(gain, rn, n_meas, _x, err=True):
	#define parameters:
	nx = _x["nx"]
	ny = _x["ny"]
	sigma_stray = _x["sigma_stray"] #e-
	sigma_dc = _x["sigma_dc"] #e-/s
	t = _x["t_rn"]
	
	#Define design variables
	#n_meas number of measurements
	
	sigma_si = gain * math.sqrt(rn**2 + (sigma_dc*t)**2 + (sigma_stray*t)**2)
	
	n = nx * ny * n_meas
	random_var = (2/n) * sigma_si**4
	random_sigma = math.sqrt(random_var)
	random = norm.rvs(scale = random_sigma)
	
	if err:
		y = sigma_si + random
	else:
		y = sigma_si
	return y
	
"""
helper fn for dark current experiment
"""
def dark_current_time_fn(i, tmin, dmax, dpow, dnum):
	i = i+1
	if i > dnum or i < 1:
		print("Bad dark current times! Sheesh!", i, dmax)
		#quit
	b = (dmax - tmin*(dnum**dpow))/(1-dnum**dpow)
	a = tmin - b
	return a * (i**dpow) + b
	
"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def dark_current_exp(gain, rn, dc, d_num, d_max, d_pow, _x, err=True):
	#define parameters:
	mu_stray = _x["mu_stray"]
	sigma_stray = _x["sigma_stray"]
	sigma_dc = _x["sigma_dc"]
	t_0 = _x["t_0"]
	nx = _x["nx"]
	ny = _x["ny"]
	
	#define design variables
	#d_num total number of exposures
	#d_max highest exposure time in the list
	#d_pow exponent of the power function describing exp times
	
	#handle design-derived variables
	t_list = []
	for i in range(d_num):
		t = dark_current_time_fn(i, tmin=t_0, dmax=d_max, dpow=d_pow, dnum=d_num)
		t_list.append(t)
		
	#Calculate mean signal at each datapoint
	signal_list = [0]*len(t_list)
	for i,_ in enumerate(signal_list):
		signal_list[i] = gain * (dc + mu_stray) * t_list[i]
	
	#Apply error to each datapoint - you could do it at the end, but why not here.
	if err:
		for i,_ in enumerate(signal_list):
			ti = t_list[i]
			random_var = gain**2 * (rn**2 + (sigma_dc*ti)**2 + (sigma_stray*ti)**2) / np.sqrt(nx**2 * ny**2)
			random_sigma = math.sqrt(random_var)
			signal_list[i] += norm.rvs(scale = random_sigma) #apply error to mean
		
	#Calculate slope (assume an intercept but don't worry abot it)
	Sxy, Sx2 = 0,0
	xmean = np.mean(t_list)
	ymean = np.mean(signal_list)
	for i,xi in enumerate(t_list):
		yi = signal_list[i]
		Sxy += (xi - xmean)*(yi - ymean)
		Sx2 += (xi - xmean)**2
	
	m = Sxy / Sx2
	return m
	
	
#This experiment does not model "photon losses due to the gate structures, electron 
#recombination within the bulk silicon itself, surface reflection, and, for very long or 
#short wavelengths, losses due to the almost complete lack of absorption by the CCD"
#- Howell, Handbook of CCD Astronomy - instead, it models how we might measure intrinsic QE
def quantum_efficiency_exp(qe, gain, rn, n_qe, t_qe, I_qe, _x, err=True):
	#define parameters
	S_pd = _x["S_pd"] #functional
	S_pd_err = _x["S_pd_err"]
	sigma_dc = _x["sigma_dc"]
	h = 6.62607015e-34 #J*Hz−1 #Planck's constant
	c = 299792458 #m/s #speed of light
	
	#define design variables
	#n_qe number of evenly-separated measurement points
	#t_qe exposure time
	#I_qe photocurrent of light source
	
	#Derived values
	if err:
		S_pd_sample = noise_to_functional(S_pd, S_pd_err)
	else:
		#need a trick for getting the mean of the functional
		0
	measure_pts = np.linspace(S_pd.xmin, S_pd.xmax, n_qe)
	
	Signal_measure = []
	for lambda_i in measure_pts:
		#see Krishnamurthy et al. 2017
		power_pd = I_qe / S_pd_sample.f(lambda_i)
		energy = h*c / (lambda_i*1e-9)
		photon_rate = power_pd / energy	
		num_photons = photon_rate * t_qe
		Signal = gain * qe.f(lambda_i) * num_photons
		
		Signal_measure.append(Signal)

	measurement = Functional(measure_pts, Signal_measure)
	measurement.spline_interp(3)
	measurement.set_xlim(qe.xmin, qe.xmax)
	measurement.set_ylim(0, measurement.ymax*1.5)
	
	if err:
		error = math.sqrt(rn**2 + (sigma_dc*t_qe)**2)
		meas_with_error = noise_to_functional(measurement, error)
		return meas_with_error
	else:
		return measurement
