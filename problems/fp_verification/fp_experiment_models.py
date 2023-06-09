import math
import numpy as np
from scipy.stats import norm, poisson

"""
Full matrix experiment model
"""
def fp_likelihood_fn(theta, d, x):
	#define interest params:
	gain = theta["gain"] #change this?
	rn = theta["rn"]
	dc = theta["dc"]
	#define design variables:
	t_gain = d["t_gain"]
	I_gain = d["I_gain"]
	n_meas_rn = d["n_meas_rn"]
	d_num = d["d_num"]
	d_max = d["d_max"]
	d_pow = d["d_pow"]
	#just pass along entire x
	
	y1 = gain_exp(gain, rn, dc, t_gain, I_gain, x)
	y2 = read_noise_exp(gain, rn, n_mean_rn, x)
	y3 = dark_current_exp(gain, rn, dc, d_num, d_max, d_pow, x)
	
	return [y1, y2, y3]


"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def gain_exp(gain, rn, dc, t, I, _x):
	#define parameters:
	P_signal = _x["P_signal"] #probability of a signal event being correctly identified as an event; 1-P_signal is false negative rate
	P_noise = _x["P_noise"] #probability of noise being incorrectly identified as an event; P_noise is false positive rate
	sigma_dc = _x["sigma_dc"]
	T_ccd = _x["T_ccd"]
	E0 = _x["E0"]
	sigma_E = _x["sigma_E"]
	w = _x["w"]
	rate_cd109 = _x["rate_cd109"] #look at cd-109 sample
	nx = _x["nx"]
	ny = _x["ny"]
	grade_size = _x["grade_size"] #we are considering 3x3 event grade sizes, based on the physics of the experiment
	
	#define design variables
	#t length of sample exposure
	#I number of exposures
	
	#Sample from Poisson distribution of # particle events
	lambda_emis = rate_cd109 * t
	N_t = poisson.rvs(lambda_emis)
	
	#calculate derived values
	grade_area = grade_size**2
	n_signal = N_t * P_signal
	n_noise = (nx*ny/grade_area) * P_noise
	n = n_signal + n_noise
	p = n_signal / n
	
	#calculate mixture distribution parameters, X = pS + (1-p)N
	mu_x = p*E0*gain/w
	var_x = p * (sigma_E/w)**2 + rn**2 + sigma_dc**2 + p*(1-p)*(E0/w)**2
	sigma_x = math.sqrt(var_x)
	
	random = norm.rvs(scale = sigma_x)
	y = mu_x + random
	return y
	

"""
theta: [0] gain [1] read noise [2] dark current
d: [0] number exposures
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def read_noise_exp(gain, rn, n_meas, _x):
	#define parameters:
	nx = _x["nx"]
	ny = _x["ny"]
	sigma_stray = _x["sigma_stray"] #e-
	sigma_dc = _x["sigma_dc"] #e-/s
	t = _x["t_rn"]
	
	#Define design variables
	#n_meas number of measurements
	
	sigma_si = gain * math.sqrt(rn**2 + (sigma_dc*t)**2 + sigma_stray**2)
	
	n = nx * ny * n_measurements
	random_var = (2/n) * sigma_si**4
	random_sigma = math.sqrt(random_var)
	random = norm.rvs(scale = random_sigma)
	
	y = sigma_si + random
	return y
	
	
"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def dark_current_exp(gain, rn, dc, _d, _x):
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

	for x in range(d_num):
		C = (d_max - math.exp(d_pow)) / (d_num-1)
		t = C * x**d_pow
		t_list.append(t)
	t_list[0] = t_0 #clobber; 100ms baseline exposure assumed
	
	t_2sum, t_1halfsum, t_3halfsum = 0
	for t in t_list:
		t_2sum += t**2
		t_1halfsum += t**(.5)
		t_3halfsum += t**(1.5)
		
	#Start calculating slope
	dc_sum, stray_sum = 0
	for t in t_list:
		dc_sum += dc * t**2
		stray_sum += mu_stray * t
		
	m = dc_sum/t_2sum + stray_sum/t_2sum

	#calc error
	random_var = (rn**2 * t_1halfsum + sigma_dc**2*t_3halfsum + sigma_stray**2 * t_1halfsum) / (nx**2 * ny**2 * t_2sum) #variance, not stddev
	random_sigma = math.sqrt(random_var)
	random = norm.rvs(scale = random_sigma)

	y = m + random
	return y