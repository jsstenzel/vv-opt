#import argparse
import os
import sys
import shutil
import csv
import fileinput

import matplotlib.pyplot as plt
from astropy.io import fits
import scipy

import numpy as np
import itertools
import multiprocessing as mp
import math
from copy import deepcopy
import scipy.optimize as optimization

os.environ['COATINGS_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/"
os.environ['TEST_DATA_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/TEST DATA/"

sys.path.insert(0, "..")
#from problems.llamas_snr import *
from approx.regression_models import *

__max_lambda_pts = 1000 #There is some loss here, but thats acceptable

"""
Full matrix experiment model
"""
def snr_likelihood_fn(theta, d, x, prior_mean, err=True):
	#define interest params:
	if type(theta) is not dict:
		print("oh howd that happen now?")
		sys.exit()
	#define design variables, gain rn dc:
	t_gain = d["t_gain"]
	I_gain = d["I_gain"]
	n_meas_rn = d["n_meas_rn"]
	d_num = d["d_num"]
	d_max = d["d_max"]
	d_pow = d["d_pow"]
	#quantum efficiency:
	n_qe = d["n_qe"]
	t_qe = d["t_qe"]

	#_wave_max = 975.0 red end
	#_wave_redgreen = 690.0
	#_wave_greenblue = 480.0
	#_wave_min = 350.0 blue end
	red_max = x["wave_max"]
	red_min = x["wave_redgreen"]
	gre_max = x["wave_redgreen"]
	gre_min = x["wave_greenblue"]
	blu_max = x["wave_greenblue"]
	blu_min = x["wave_min"]
	
	y_gain_red = gain_exp(theta["gain_red"], theta["rn_red"], theta["dc_red"], t_gain, I_gain, x, prior_mean["gain_red"], err)
	y_gain_gre = gain_exp(theta["gain_gre"], theta["rn_gre"], theta["dc_gre"], t_gain, I_gain, x, prior_mean["gain_gre"], err)
	y_gain_blu = gain_exp(theta["gain_blu"], theta["rn_blu"], theta["dc_blu"], t_gain, I_gain, x, prior_mean["gain_blu"], err)
	y_gain = [y_gain_red, y_gain_gre, y_gain_blu]
	
	y_rn_red = read_noise_exp(theta["gain_red"], theta["rn_red"], n_meas_rn, x, prior_mean["rn_red"], err)
	y_rn_gre = read_noise_exp(theta["gain_gre"], theta["rn_gre"], n_meas_rn, x, prior_mean["rn_gre"], err)
	y_rn_blu = read_noise_exp(theta["gain_blu"], theta["rn_blu"], n_meas_rn, x, prior_mean["rn_blu"], err)
	y_rn = [y_rn_red, y_rn_gre, y_rn_blu]
	
	y_dc_red = dark_current_exp(theta["gain_red"], theta["rn_red"],  theta["dc_red"], d_num, d_max, d_pow, x, prior_mean["dc_red"], err)
	y_dc_gre = dark_current_exp(theta["gain_gre"], theta["rn_gre"],  theta["dc_gre"], d_num, d_max, d_pow, x, prior_mean["dc_gre"], err)
	y_dc_blu = dark_current_exp(theta["gain_blu"], theta["rn_blu"],  theta["dc_blu"], d_num, d_max, d_pow, x, prior_mean["dc_blu"], err)
	y_dc = [y_dc_red, y_dc_gre, y_dc_blu]
	
	y_qe_red_t = quantum_efficiency_exp(theta["qe_red_t"], theta["gain_red"], theta["rn_red"], n_qe, t_qe, red_min, red_max, red_max-blu_min, x, prior_mean["qe_red_t"], err)
	y_qe_gre_t = quantum_efficiency_exp(theta["qe_gre_t"], theta["gain_gre"], theta["rn_gre"], n_qe, t_qe, gre_min, gre_max, red_max-blu_min, x, prior_mean["qe_gre_t"], err)
	y_qe_blu_t = quantum_efficiency_exp(theta["qe_blu_t"], theta["gain_blu"], theta["rn_blu"], n_qe, t_qe, blu_min, blu_max, red_max-blu_min, x, prior_mean["qe_blu_t"], err)
	
	y_vph_red_t = measure_thru_vph(theta["vph_red_t"], d["d_vph_n_pts"], red_min, red_max, x["vph_meas_stddev"], prior_mean["vph_red_t"], err)
	y_vph_gre_t = measure_thru_vph(theta["vph_gre_t"], d["d_vph_n_pts"], gre_min, gre_max, x["vph_meas_stddev"], prior_mean["vph_gre_t"], err)
	y_vph_blu_t = measure_thru_vph(theta["vph_blu_t"], d["d_vph_n_pts"], blu_min, blu_max, x["vph_meas_stddev"], prior_mean["vph_blu_t"], err)
	
	y_sl_t = measure_thru_sigmoid(theta["sl_t"], d["d_dichroic_n_pts"], x["wave_min"], x["wave_max"], x["sl_meas_stddev"], prior_mean["sl_t"], err)
	y_bg_t = measure_thru_sigmoid(theta["bg_t"], d["d_dichroic_n_pts"], x["wave_min"], x["wave_max"], x["bg_meas_stddev"], prior_mean["bg_t"], err)
	
	#collimator
	y_coll_t = thru_measurement([theta["coll_t"]], d["d_coll_n_pts"], x["wave_min"], x["wave_max"], x["coll_meas_stddev"], [prior_mean["coll_t"]], err)
	
	#lenses
	#TODO get rid of prism, replace this with 3 experiments that each measure the total thru of a full camera
	red_lenses = [theta["red_l1_t"],theta["red_l2_t"],theta["red_l3_t"],theta["red_l4_t"],theta["red_l5_t"],theta["red_l6_t"],theta["red_l7_t"]]
	red_priormeans = [prior_mean["red_l1_t"],prior_mean["red_l2_t"],prior_mean["red_l3_t"],prior_mean["red_l4_t"],prior_mean["red_l5_t"],prior_mean["red_l6_t"],prior_mean["red_l7_t"]]
	y_red_cam = thru_measurement(red_lenses, d["d_redcam_n_pts"], red_min, red_max, x["lens_meas_err"], red_priormeans, err)
	
	gre_lenses = [theta["gre_l1_t"],theta["gre_l2_t"],theta["gre_l3_t"],theta["gre_l4_t"],theta["gre_l5_t"],theta["gre_l6_t"],theta["gre_l7_t"]]
	gre_priormeans = [prior_mean["gre_l1_t"],prior_mean["gre_l2_t"],prior_mean["gre_l3_t"],prior_mean["gre_l4_t"],prior_mean["gre_l5_t"],prior_mean["gre_l6_t"],prior_mean["gre_l7_t"]]
	y_gre_cam = thru_measurement(gre_lenses, d["d_greencam_n_pts"], gre_min, gre_max, x["lens_meas_err"], gre_priormeans, err)
	
	blu_lenses = [theta["blu_l1_t"],theta["blu_l2_t"],theta["blu_l3_t"],theta["blu_l4_t"],theta["blu_l5_t"],theta["blu_l6_t"],theta["blu_l7_t"],theta["blu_l8_t"]]
	blu_priormeans = [prior_mean["blu_l1_t"],prior_mean["blu_l2_t"],prior_mean["blu_l3_t"],prior_mean["blu_l4_t"],prior_mean["blu_l5_t"],prior_mean["blu_l6_t"],prior_mean["blu_l7_t"],prior_mean["blu_l8_t"]]
	y_blu_cam = thru_measurement(blu_lenses, d["d_bluecam_n_pts"], blu_min, blu_max, x["lens_meas_err"], blu_priormeans, err)
	
	y_lenses = [
		y_red_cam, y_gre_cam, y_blu_cam
	]
	
	y_frd = simple_measurement(theta["fiber_frd"], x["frd_meas_err"], d["d_frd_n_meas"], prior_mean["fiber_frd"], err)
	
	y = [*y_gain, *y_rn, *y_dc, *y_qe_red_t, *y_qe_gre_t, *y_qe_blu_t, *y_vph_red_t, *y_vph_gre_t, *y_vph_blu_t, *y_sl_t, *y_bg_t, y_coll_t, *y_lenses, y_frd]
	return y

	
"""
Add standard error to a direct measurement of the input
theta: float
prior_mean: float
"""
def simple_measurement(theta, stddev_meas, n_meas, prior_mean, err=True):
	if n_meas > 1:
		y = theta
		stddev = stddev_meas / np.sqrt(n_meas)

		if err:
			random = scipy.stats.norm.rvs(scale = stddev)
			y += random

		return y
		
	else: 
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		
		return prior_mean

"""
measure the mean throughput of one or more optical elements
elem_list: list of Gaussian Process objects
prior_mean_list: list of Gaussian Process objects
"""
def thru_measurement(elem_list, n_meas, wave_min, wave_max, meas_stddev, prior_mean_list, err=True):
	if err:
		stddev = meas_stddev
	else:
		stddev = 0
		
	if n_meas > 1:
		#choose the measurement points; this +2 -2 thing is to make sure im not measuring the endpoints
		measurement_pts = np.linspace(wave_min, wave_max, num=n_meas+2)
		measurement_pts = measurement_pts[1:-1]
		
		#Measure each element, and stack all the curves on top of each other
		y_thru = []
		for elem in elem_list:
			thru_i = elem.measure(measurement_pts, stddev)
			if y_thru == []:
				y_thru = thru_i
			else:
				y_thru = [y*i for y,i in zip(y_thru, thru_i)]
	else: 
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		max_pts = np.linspace(wave_min, wave_max, num=__max_lambda_pts)
		y_thru = []
		for elem in prior_mean_list:
			thru_i = elem.measure(max_pts, 0)
			if y_thru == []:
				y_thru = thru_i
			else:
				y_thru = [y*i for y,i in zip(y_thru, thru_i)]

	return np.mean(y_thru)
	
"""
measure and fit a sigmoid to a throughput
theta_gp: Gaussian Process object
prior_mean: Gaussian Process object
"""
def measure_thru_sigmoid(theta_gp, d_meas_pts, wave_min, wave_max, meas_stddev, prior_mean, err=True):
	if err:
		stddev = meas_stddev
	else:
		stddev = 0
		
	#Handle measurement number: you can do 0, or 4 or more
	if d_meas_pts == 0:
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		max_pts = np.linspace(wave_min, wave_max, num=__max_lambda_pts)
		y_thru = prior_mean.measure(max_pts, 0)
		lval, step_pt, rval, power, _ = sigmoid_fit_throughput(max_pts, y_thru, doPlot=False, doErr=False)
		return [lval, step_pt, rval, power]
	elif d_meas_pts < 4:
		num_pts = 4 #I could do smarter stuff here
	else:
		num_pts = d_meas_pts
	
	#choose the measurement points
	measurement_pts = np.linspace(wave_min, wave_max, num=num_pts)
	#print(measurement_pts, flush=True)
	
	#make the measurements, assuming that there is one y_i for each measurement point ki
	y_thru = theta_gp.measure(measurement_pts, stddev)
	#print(y_thru, flush=True)
	
	#apply the 0..1 boundaries
	for i,yi in enumerate(y_thru):
		if yi < 0:
			y_thru[i] = 0
		if yi > 1:
			y_thru[i] = 1
			
	#convert measurements back to y params
	lval, step_pt, rval, power, _ = sigmoid_fit_throughput(measurement_pts, y_thru, doPlot=False, doErr=False)
	
	return [lval, step_pt, rval, power]
	
"""
Measure and fit a parabola to a throughput
theta_gp: Gaussian Process object
prior_mean: Gaussian Process object
"""
def measure_thru_vph(theta_gp, d_meas_pts, wave_min, wave_max, meas_stddev, prior_mean, err=True):
	if err:
		stddev = meas_stddev
	else:
		stddev = 0
		
	#Handle measurement number: you can do 0, or 3 or more
	if d_meas_pts == 0:
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		max_pts = np.linspace(wave_min, wave_max, num=__max_lambda_pts)
		y_thru = prior_mean.measure(max_pts, 0)
		popt, _ = poly_fit_throughput(max_pts, y_thru, 2, doPlot=False, doErr=False, doCov=False)
		return popt
	elif d_meas_pts < 3:
		num_pts = 3 #I could do smarter stuff here
	else:
		num_pts = d_meas_pts
	
	#choose the measurement points; this +2 -2 thing is to make sure im not measuring the endpoints
	measurement_pts = np.linspace(wave_min, wave_max, num=num_pts+2)
	measurement_pts = measurement_pts[1:-1]
	
	#make the measurements
	y_thru = theta_gp.measure(measurement_pts, stddev)
	
	#apply the 0..1 boundaries
	for i,yi in enumerate(y_thru):
		if yi < 0:
			y_thru[i] = 0
		if yi > 1:
			y_thru[i] = 1
			
	#convert measurements back to y params
	#NOTE: we assume, historically, that the VPH throughput is parabolic.
	popt, _ = poly_fit_throughput(measurement_pts, y_thru, 2, doPlot=False, doErr=False, doCov=False)
		
	return popt
	
	
"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
Assumptions:
* All decays of the sample make it to the CCD
* Signal and noise events will never overlap
* w has no uncertainty
"""
def gain_exp(gain, rn, dc, t, I, _x, prior_mean, err=True):	
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
	
	if I <= 0 or t <= 0:
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		p = 1 #this is what happens if you assume no false positives, perfect experiment, just get the right units
		mu_x = p*E0*prior_mean/w	
		return mu_x
	
	#calculate derived values
	grade_area = grade_size**2
	n_signal = N_decays * P_signal
	n_noise = (nx*ny/grade_area) * P_noise
	n = n_signal + n_noise
	p = n_signal / n
	
	#calculate mixture distribution parameters, X = pS + (1-p)N
	mu_x = p*E0*gain/w	
	
	if err:
		var_x = p * (sigma_E/w)**2 + rn**2 + sigma_dc**2 + p*(1-p)*(E0/w)**2
		sigma_x = math.sqrt(var_x)
		stddev = sigma_x/np.sqrt(I*n)
		random = scipy.stats.norm.rvs(scale = stddev)
		y = mu_x + random
	else:
		y = mu_x
	return y
	

"""
theta: [0] gain [1] read noise [2] dark current
d: [0] number exposures
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def read_noise_exp(gain, rn, n_meas, _x, prior_mean, err=True):
	if n_meas <= 0:
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		return gain * prior_mean
		
	#define parameters:
	nx = _x["nx"]
	ny = _x["ny"]
	sigma_stray = _x["sigma_stray"] #e-
	sigma_dc = _x["sigma_dc"] #e-/s
	t = _x["t_rn"]
	
	#Define design variables
	#n_meas number of measurements
	
	sigma_si = gain * math.sqrt(rn**2 + (sigma_dc*t)**2 + (sigma_stray*t)**2)
	
	if err:
		n = nx * ny * n_meas
		random_var = (2/n) * sigma_si**4
		random_sigma = math.sqrt(random_var)
		random = scipy.stats.norm.rvs(scale = random_sigma)
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
def dark_current_exp(gain, rn, dc, d_num, d_max, d_pow, _x, prior_mean, err=True):
	if d_num <= 1:
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		return prior_mean
		
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
			signal_list[i] += scipy.stats.norm.rvs(scale = random_sigma) #apply error to mean
		
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
	
"""
Measure the signal of a QE experiment, and fit a linear regression with fourier elements to the signal
qe: Gaussian Process object
prior_mean: Gaussian Process object
"""
def quantum_efficiency_exp(theta_qe, gain, rn, n_qe, t_qe, wave_min, wave_max, full_bandpass, _x, prior_mean, err=True):				
	#define parameters
	S_pd = _x["S_pd"] #functional
	S_pd_err = _x["S_pd_meas_err"]
	sigma_dc = _x["sigma_dc"]
	spectral_power = _x["spectral_power"] #W / nm
	h = 6.62607015e-34 #J*Hz−1 #Planck's constant
	c = 299792458 #m/s #speed of light
	
	#Handle measurement number: you can do 0, or 3 or more
	if n_qe <= 0:
		#This means we don't run the experiment
		#Model the imputation that would occur: 
		# - Ignore the provided theta and assume theta is the mean of the prior
		# - calculate the y that would result in
		qe = prior_mean
		num_pts = __max_lambda_pts
		err = False
	elif n_qe < 3:
		num_pts = 3
		qe = theta_qe
	else:
		num_pts = n_qe
		qe = theta_qe
	
	#define design variables
	#n_qe number of evenly-separated measurement points
	#t_qe exposure time
	#I_qe photocurrent of light source - spectral power * wavelength, see Krishnamurthy et al. 2016
	
	#times smaller than this just aren't meaningful
	if t_qe <= 0.1:
		t = 0.1
	else:
		t = t_qe
	
	#Take sample measurements of source and qe
	measure_pts = np.linspace(wave_min, wave_max, num_pts)
	
	if err:
		S_pd_sample = S_pd.measure(measure_pts, S_pd_err)
	else:
		S_pd_sample = S_pd.measure(measure_pts, 0)

	qe_sample = qe.measure(measure_pts, 0)
	
	Signal_measure = []
	for lambda_i, qe_i, S_i in zip(measure_pts, qe_sample, S_pd_sample):
		#see Krishnamurthy et al. 2017
		I_qe = spectral_power * lambda_i
		power_pd = I_qe / S_i
		energy = h*c / (lambda_i*1e-9)
		photon_rate = power_pd / energy	
		num_photons = photon_rate * t
		Signal = gain * qe_i * num_photons
		
		if err:
			ccd_error = math.sqrt(rn**2 + (sigma_dc*t)**2)
			Signal += scipy.stats.norm.rvs(scale = ccd_error)
		
		Signal_measure.append(Signal)
		
	#convert measurements back to y params
	coeffs, inter, _, _ = linreg_fourier_throughput(measure_pts, Signal_measure, 2, full_bandpass, doPlot=False, doErr=False)
		
	return [inter, coeffs[0], coeffs[1], coeffs[2], coeffs[3]]