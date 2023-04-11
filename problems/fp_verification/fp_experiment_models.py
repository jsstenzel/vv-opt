import math
import numpy as np
from scipy.stats import norm, poisson

"""
Full matrix experiment model
"""
def fp_likelihood_fn(theta, d):
	#define interest params:
	#gain = theta[0]
	#rn = theta[1]
	#dc = theta[2]
    #define design vars:
    
    y1 = read_noise_exp([theta[0],theta[1]], d)
    y2 = dark_current_exp([theta[0],theta[1],theta[2]], d)
    y3 = gain_exp([theta[0],theta[1],theta[2]], d)
    
    return [y1, y2, y3]


"""
theta: [0] gain [1] read noise [2] dark current
d: [0] number exposures
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def read_noise_exp(_theta, _d):
	#define interest params:
	gain = _theta[0]
	rn = _theta[1]
	#define parameters:
	nx = 2048
	ny = 2048
	sigma_stray = .1 #e-
	sigma_dc = .5 #e-/s
	t = .1 #100 ms exposure, best for both noise and 
	
	sigma_si = gain * math.sqrt(rn**2 + (sigma_dc*t)**2 + sigma_stray**2)
	
	n = nx * ny * _d[0]
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
def dark_current_exp(_theta, _d):
	#define interest params:
	gain = _theta[0]
	rn = _theta[1]
	dc = _theta[2]
	#define parameters:
	mu_stray = 0
	sigma_stray = .1 #e-
	sigma_dc = .5 #e-/s
	t_0 = 0.1 #100ms baseline exposure assumed
	nx = 2048
	ny = 2048
	#define design variables
	d_num = _d[0] #total number of exposures
	d_max = _d[1] #highest exposure time in the list
	d_pow = _d[2] #exponent of the power function describing exp times
	
	#handle design-derived variables
	t_list = [0]*d_num

	for x,t in enumerate(t_list):
		C = (d_max - math.exp(d_pow)) / (d_num-1)
		t = C * x**d_pow
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
	
	
"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def gain_exp(_theta, _d):
	#define interest params:
	gain = _theta[0]
	rn = _theta[1]
	dc = _theta[2]
	#define parameters:
	P_signal = #probability of a signal event being correctly identified as an event; 1-P_signal is false negative rate
	P_noise = #probability of noise being incorrectly identified as an event; P_noise is false positive rate
	sigma_dc = .5 #e-/s
	T_ccd =
	E0 =
	sigma_E =
	w =
	rate_cd109 = #look at cd-109 sample
	nx = 2048
	ny = 2048
	grade_size = 3 #we are considering 3x3 event grade sizes, based on the physics of the experiment
	#define design variables
	t = _d[0]
	I = _d[1] #number of iterations
	
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
	