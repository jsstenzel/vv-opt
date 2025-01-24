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
import more_itertools
import multiprocessing as mp
import math
from copy import deepcopy
import scipy.optimize as optimization

os.environ['COATINGS_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/"
os.environ['TEST_DATA_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/TEST DATA/"

sys.path.insert(0, "..")
#from problems.llamas_snr import *
from approx.regression_models import *


"""
Full matrix experiment model
"""
def snr_likelihood_fn(theta, d, x, err=True):
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
	
	y_gain_red = gain_exp(theta["gain_red"], theta["rn_red"], theta["dc_red"], t_gain, I_gain, x, err)
	y_gain_gre = gain_exp(theta["gain_gre"], theta["rn_gre"], theta["dc_gre"], t_gain, I_gain, x, err)
	y_gain_blu = gain_exp(theta["gain_blu"], theta["rn_blu"], theta["dc_blu"], t_gain, I_gain, x, err)
	y_gain = [y_gain_red, y_gain_gre, y_gain_blu]
	
	y_rn_red = read_noise_exp(theta["gain_red"], theta["rn_red"], n_meas_rn, x, err)
	y_rn_gre = read_noise_exp(theta["gain_gre"], theta["rn_gre"], n_meas_rn, x, err)
	y_rn_blu = read_noise_exp(theta["gain_blu"], theta["rn_blu"], n_meas_rn, x, err)
	y_rn = [y_rn_red, y_rn_gre, y_rn_blu]
	
	y_dc_red = dark_current_exp(theta["gain_red"], theta["rn_red"],  theta["dc_red"], d_num, d_max, d_pow, x, err)
	y_dc_gre = dark_current_exp(theta["gain_gre"], theta["rn_gre"],  theta["dc_gre"], d_num, d_max, d_pow, x, err)
	y_dc_blu = dark_current_exp(theta["gain_blu"], theta["rn_blu"],  theta["dc_blu"], d_num, d_max, d_pow, x, err)
	y_dc = [y_dc_red, y_dc_gre, y_dc_blu]
	
	qe_red_params = [theta["qe_red_t0"], theta["qe_red_t1"], theta["qe_red_t2"], theta["qe_red_t3"], theta["qe_red_t4"]]
	qe_gre_params = [theta["qe_gre_t0"], theta["qe_gre_t1"], theta["qe_gre_t2"], theta["qe_gre_t3"], theta["qe_gre_t4"]]
	qe_blu_params = [theta["qe_blu_t0"], theta["qe_blu_t1"], theta["qe_blu_t2"], theta["qe_blu_t3"], theta["qe_blu_t4"]]
	y_qe_red_t = quantum_efficiency_exp(qe_red_params, theta["gain_red"], theta["rn_red"], n_qe, t_qe, red_min, red_max, red_max-blu_min, x, err)
	y_qe_gre_t = quantum_efficiency_exp(qe_gre_params, theta["gain_gre"], theta["rn_gre"], n_qe, t_qe, gre_min, gre_max, red_max-blu_min, x, err)
	y_qe_blu_t = quantum_efficiency_exp(qe_blu_params, theta["gain_blu"], theta["rn_blu"], n_qe, t_qe, blu_min, blu_max, red_max-blu_min, x, err)
	
	vph_red_params = [theta["vph_red_t0"], theta["vph_red_t1"], theta["vph_red_t2"]]#, theta["vph_red_t3"]]
	vph_gre_params = [theta["vph_gre_t0"], theta["vph_gre_t1"], theta["vph_gre_t2"]]#, theta["vph_gre_t3"]]
	vph_blu_params = [theta["vph_blu_t0"], theta["vph_blu_t1"], theta["vph_blu_t2"]]#, theta["vph_blu_t3"]]
	y_vph_red_t = measure_thru_vph(vph_red_params, d["d_vph_n_pts"], red_min, red_max, x["vph_meas_stddev"], err)
	y_vph_gre_t = measure_thru_vph(vph_gre_params, d["d_vph_n_pts"], gre_min, gre_max, x["vph_meas_stddev"], err)
	y_vph_blu_t = measure_thru_vph(vph_blu_params, d["d_vph_n_pts"], blu_min, blu_max, x["vph_meas_stddev"], err)
	
	sl_params = [theta["sl_t0"], theta["sl_t1"], theta["sl_t2"], theta["sl_t3"]]
	bg_params = [theta["bg_t0"], theta["bg_t1"], theta["bg_t2"], theta["bg_t3"]]
	y_sl_t = measure_thru_sigmoid(sl_params, d["d_dichroic_n_pts"], x["wave_min"], x["wave_max"], x["sl_meas_stddev"], err)
	y_bg_t = measure_thru_sigmoid(bg_params, d["d_dichroic_n_pts"], x["wave_min"], x["wave_max"], x["bg_meas_stddev"], err)
	
	y_frd = simple_measurement(theta["fiber_frd"], x["frd_meas_err"], d["d_frd_n_meas"], err)
	
	y = [*y_gain, *y_rn, *y_dc, *y_qe_red_t, *y_qe_gre_t, *y_qe_blu_t, *y_vph_red_t, *y_vph_gre_t, *y_vph_blu_t, *y_sl_t, *y_bg_t, y_frd]
	return y

	
#Adding standard error to a direct measurement of the input
def simple_measurement(theta, stddev_meas, n_meas, err=True):
	y = theta
	stddev = stddev_meas / np.sqrt(n_meas)

	if err:
		random = scipy.stats.norm.rvs(scale = stddev)
		y += random

	return y

	
def measure_thru_sigmoid(params, d_meas_pts, wave_min, wave_max, meas_stddev, err=True):
	if err:
		stddev = meas_stddev
	else:
		stddev = 0
		
	#convert the theta params to a GP
	lambda_pts = np.linspace(wave_min, wave_max, num=math.ceil((wave_max-wave_min)/0.1))
	#print(lambda_pts, flush=True)
	thru_pts = throughput_from_sigmoidfit_coeffs(params[0], params[1], params[2], params[3], lambda_pts)
	#print(thru_pts, flush=True)
	gp_throughput = define_functional(lambda_pts, thru_pts, order=1)
	
	#choose the measurement points
	measurement_pts = np.linspace(wave_min, wave_max, num=d_meas_pts)
	#print(measurement_pts, flush=True)
	
	#make the measurements, assuming that there is one y_i for each measurement point ki
	y_thru = gp_throughput.eval_gp_cond(measurement_pts, stddev)
	#print(y_thru, flush=True)
	
	#apply the 0..1 boundaries
	for i,yi in enumerate(y_thru):
		if yi < 0:
			y_thru[i] = 0
		if yi > 1:
			y_thru[i] = 1
			
	#convert measurements back to y params
	lval, step_pt, rval, power = sigmoid_fit_throughput(measurement_pts, y_thru, doPlot=False, doErr=False)
	
	return [lval, step_pt, rval, power]
	
	
def measure_thru_vph(params, d_meas_pts, wave_min, wave_max, meas_stddev, err=True):
	if err:
		stddev = x["vph_meas_stddev"]
	else:
		stddev = 0
		
	#convert the theta params to a GP
	lambda_pts = np.linspace(wave_min, wave_max, num=math.ceil((wave_max-wave_min)/0.1))
	thru_pts = throughput_from_polyfit_coeffs(params, lambda_pts)
	gp_throughput = define_functional(lambda_pts, thru_pts, order=1)
	
	#choose the measurement points
	measurement_pts = np.linspace(wave_min, wave_max, num=d_meas_pts+2)
	measurement_pts = measurement_pts[1:-1]
	
	#make the measurements
	y_thru = gp_throughput.eval_gp_cond(measurement_pts, stddev)
	
	#apply the 0..1 boundaries
	for i,yi in enumerate(y_thru):
		if yi < 0:
			y_thru[i] = 0
		if yi > 1:
			y_thru[i] = 1
			
	#convert measurements back to y params
	#NOTE: we assume, historically, that the VPH throughput is parabolic.
	popt = poly_fit_throughput(measurement_pts, y_thru, 2, doPlot=False, doErr=False)
		
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
	stddev = sigma_x/np.sqrt(math.sqrt(I)*n)
	
	random = scipy.stats.norm.rvs(scale = stddev)
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
	random = scipy.stats.norm.rvs(scale = random_sigma)
	
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
	
	
#This experiment does not model "photon losses due to the gate structures, electron 
#recombination within the bulk silicon itself, surface reflection, and, for very long or 
#short wavelengths, losses due to the almost complete lack of absorption by the CCD"
#- Howell, Handbook of CCD Astronomy - instead, it models how we might measure intrinsic QE
def quantum_efficiency_exp(params, gain, rn, n_qe, t_qe, wave_min, wave_max, full_bandpass, _x, err=True):
	#define parameters
	S_pd = _x["S_pd"] #functional
	S_pd_err = _x["S_pd_meas_err"]
	sigma_dc = _x["sigma_dc"]
	spectral_power = _x["spectral_power"] #W / nm
	h = 6.62607015e-34 #J*Hz−1 #Planck's constant
	c = 299792458 #m/s #speed of light
	
	#choose the measurement points
	measure_pts = np.linspace(wave_min, wave_max, n_qe)
	
	#get the qe at each measurement point
	qe_sample = throughput_from_linfourier_coeffs(params[1:], params[0], 2, wave_max-wave_min, measure_pts)
	
	#apply the 0..1 boundaries to the qe
	for i,yi in enumerate(qe_sample):
		if yi < 0:
			qe_sample[i] = 0
		if yi > 1:
			qe_sample[i] = 1
	
	if err:
		S_pd_sample = S_pd.eval_gp_cond(measure_pts, S_pd_err)
	else:
		S_pd_sample = S_pd.eval_gp_cond(measure_pts, 0)
	
	Signal_measure = []
	for lambda_i, qe_i, S_i in zip(measure_pts, qe_sample, S_pd_sample):
		#see Krishnamurthy et al. 2017
		I_qe = spectral_power * lambda_i
		power_pd = I_qe / S_i #I could add error here, associated with the error of the photodiode
		
		energy = h*c / (lambda_i*1e-9)
		photon_rate = power_pd / energy	
		num_photons = photon_rate * t_qe
		Signal = (qe_i / gain) * num_photons
		
		if err:
			ccd_error = math.sqrt(rn**2 + (sigma_dc*t_qe))
			Signal += scipy.stats.norm.rvs(scale = ccd_error)
		
		Signal_measure.append(Signal)
		
	
	#convert measurements back to y params
	coeffs, inter, _, _ = linreg_fourier_throughput(measure_pts, Signal_measure, 2, full_bandpass, doPlot=False, doErr=False)
		
	return [inter, coeffs[0], coeffs[1], coeffs[2], coeffs[3]]
	
	
	
def cost_model(theta, x):
	return 0