#import argparse
import os
import sys
import shutil
import csv
import fileinput

import numpy as np
import matplotlib.pyplot as plt
import itertools
import multiprocessing as mp
import math
from copy import deepcopy

sys.path.insert(0, "..")
sys.path.append('llamas-etc')

import spectrograph as spec
import observe
import llamas_plot_throughput
from approx.regression_models import *


#os.environ['COATINGS_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/"
#os.environ['TEST_DATA_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/TEST DATA/"

################################
#Utility functions
################################

def build_model(_dir):
	os.environ["COATINGS_PATH"] =_dir

	llamas_red = spec.Spectrograph('LLAMAS_RED')
	llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
	llamas_green = spec.Spectrograph('LLAMAS_GREEN') 

	llamas_red.build_model('llamas-etc/llamas_red1.def')
	llamas_blue.build_model('llamas-etc/llamas_blue1.def')
	llamas_green.build_model('llamas-etc/llamas_green1.def')
	#llamas_waves = np.array(np.concatenate([llamas_blue.waves,llamas_green.waves,llamas_red.waves]))
	llamas_waves = llamas_red.waves
	
	spectrograph_models = [llamas_red, llamas_blue, llamas_green, llamas_waves]
	return spectrograph_models
	
def thru(gp):
	thru_raw = gp.evaluate()
	thru_get = []
	for t in thru_raw:
		if t >= 1.0:
			thru_get.append(1.0)
		elif t <= 0.0:
			thru_get.append(0.0)
		else:
			thru_get.append(t)
	return thru_get


################################
#Utility functions
################################

#Wrapper function for model call
#def sensitivity_hlva(x, d, col, config_table, spectrograph_models):
def sensitivity_hlva(theta, x, verbose=True):
	spect_models = x["spect_model_tracker"]
	llamas_red = spect_models[0]
	llamas_blue = spect_models[1]
	llamas_green = spect_models[2]
	waves = spect_models[3]
	
	
	##############################
	###Break down x
	##############################
	
	#focal plane
	nx = x["nx"]
	ny = x["ny"]
	pixel_size_mm = x["pixel_size_mm"]
	#llamas
	resolution = x["resolution"]
	wave_min = x["wave_min"]
	wave_max = x["wave_max"]
	f_col = x["f_collimator"]
	f_cam = x["f_camera"]
	fratio_col = x["fratio_collimator"]
	fratio_com = x["fratio_camera"]
	blade_obscure = x["blade_obscure"]
	d_beam = x["d_beam"]
	#fiber
	fiber_theta = x["fiber_theta"]
	fiber_dcore = x["fiber_dcore"]
	elem_mla = x["microlens"]#"ECI_FusedSilica.txt")]
	elem_ar = x["fiber_ar"]#"fiber_ar.txt")]
	elem_fiber = x["fiber_internal"]#"Polymicro_FBPI_8m.txt")]
	#spectrograph
	#elem_collimator = x["collimator"]#"dielectric_mirror.txt")]
	elem_prism = x["prism"]#"ECI_FusedSilica.txt")]
	elem_window = x["sensor_window"]#"ECI_FusedSilica.txt")]
	elem_glass_red = x["sensor_glass_red"]#"llamas_internal_red.txt")]
	elem_glass_gre = x["sensor_glass_gre"]#"llamas_internal.txt")]
	elem_glass_blu = x["sensor_glass_blu"]#"llamas_internal_blue.txt")]
	#snr
	tau = x["tau"]
	skyspec = x["skyspec"]
	
	red_max = x["wave_max"]
	red_min = x["wave_redgreen"]
	gre_max = x["wave_redgreen"]
	gre_min = x["wave_greenblue"]
	blu_max = x["wave_greenblue"]
	blu_min = x["wave_min"]

	##############################
	###Break down theta
	##############################
	
	rn_red = theta["rn_red"]
	rn_gre = theta["rn_gre"]
	rn_blu = theta["rn_blu"]
	dc_red = theta["dc_red"]
	dc_gre = theta["dc_gre"]
	dc_blu = theta["dc_blu"]

	lambda_pts = np.linspace(wave_min, wave_max, num=math.ceil((wave_max-wave_min)/0.1))
	
	#update qe
	elem_qe_red = theta["qe_red_t"]#CCD42-40_dd.txt
	elem_qe_gre = theta["qe_gre_t"]#CCD42-40_green.txt
	elem_qe_blu = theta["qe_blu_t"]#CCD42-40_blue.txt
	llamas_red.sensor.qe.setThroughputTab(elem_qe_red.fine_domain, thru(elem_qe_red))
	llamas_green.sensor.qe.setThroughputTab(elem_qe_gre.fine_domain, thru(elem_qe_gre))
	llamas_blue.sensor.qe.setThroughputTab(elem_qe_blu.fine_domain, thru(elem_qe_blu))
	
	#update vph
	elem_vph_red = theta["vph_red_t"]#wasach_llamas2200_red.txt
	elem_vph_gre = theta["vph_gre_t"]#wasach_llamas2200_green.txt
	elem_vph_blu = theta["vph_blu_t"]#wasach_llamas2200_blue.txt
	llamas_red.grating.blaze.setThroughputTab(elem_vph_red.fine_domain, thru(elem_vph_red))
	llamas_green.grating.blaze.setThroughputTab(elem_vph_gre.fine_domain, thru(elem_vph_gre))
	llamas_blue.grating.blaze.setThroughputTab(elem_vph_blu.fine_domain, thru(elem_vph_blu))
	
	#update collimator throughputs
	elem_collimator = theta["coll_t"]
	llamas_red.elements[0].surfaces[0].setThroughputTab(elem_collimator.fine_domain, thru(elem_collimator))
	llamas_green.elements[0].surfaces[0].setThroughputTab(elem_collimator.fine_domain, thru(elem_collimator))
	llamas_blue.elements[0].surfaces[0].setThroughputTab(elem_collimator.fine_domain, thru(elem_collimator))
			
	#update dichroic throughputs
	elem_dichroic_sl = theta["sl_t"]#ECI_FusedSilica.txt
	elem_dichroic_bg = theta["bg_t"]#ECI_FusedSilica.txt
	llamas_red.elements[1].surfaces[0].setThroughputTab(elem_dichroic_sl.fine_domain, thru(elem_dichroic_sl))
	llamas_green.elements[1].surfaces[0].setThroughputTab(elem_dichroic_sl.fine_domain, [1.0-t for t in thru(elem_dichroic_sl)])		
	llamas_green.elements[2].surfaces[0].setThroughputTab(elem_dichroic_bg.fine_domain, elem_dichroic_bg.evaluate())
	llamas_blue.elements[1].surfaces[0].setThroughputTab(elem_dichroic_sl.fine_domain, [1.0-t for t in thru(elem_dichroic_sl)])
	llamas_blue.elements[2].surfaces[0].setThroughputTab(elem_dichroic_bg.fine_domain, [1.0-t for t in thru(elem_dichroic_bg)])
	
	fiber_frd = theta["fiber_frd"] if theta["fiber_frd"]>0 else 0
	
	##############################
	###Load in the new values to the spect_models
	##############################
	
	#update spectrograph definitions
	for llamas_channel in [llamas_red, llamas_green, llamas_blue]:
		llamas_channel.R		  =	resolution
		llamas_channel.wv_min	 =	wave_min 
		llamas_channel.wv_max	 =	wave_max
		llamas_channel.f_col	  =	f_col
		llamas_channel.f_cam	  =	f_cam
		llamas_channel.Fratio_col =	fratio_col
		llamas_channel.Fratio_cam =	fratio_com
		llamas_channel.Dbeam	  =	d_beam
		llamas_channel.blade_obscuration = blade_obscure
		
	#update fiber
	for llamas_channel in [llamas_red, llamas_green, llamas_blue]:
		llamas_channel.fiber.frd_loss = fiber_frd
		llamas_channel.fiber.theta_fib = fiber_theta
		llamas_channel.fiber.dFib = fiber_dcore
		for surf in llamas_channel.fiber.elements[0].surfaces: #mla
			surf.setThroughputTab(elem_mla.fine_domain, thru(elem_mla))
		for surf in llamas_channel.fiber.elements[1].surfaces: #ar
			surf.setThroughputTab(elem_ar.fine_domain, thru(elem_ar))
		for surf in llamas_channel.fiber.elements[2].surfaces: #fiber
			surf.setThroughputTab(elem_fiber.fine_domain, thru(elem_fiber))
	
	red_elem_list = [elem_prism, theta["red_l1_t"], theta["red_l2_t"], theta["red_l3_t"], theta["red_l4_t"], theta["red_l5_t"], theta["red_l6_t"], theta["red_l7_t"], elem_window, elem_glass_red]
	gre_elem_list = [elem_prism, theta["gre_l1_t"], theta["gre_l2_t"], theta["gre_l3_t"], theta["gre_l4_t"], theta["gre_l5_t"], theta["gre_l6_t"], theta["gre_l7_t"], elem_window, elem_glass_gre]
	blu_elem_list = [elem_prism, theta["blu_l1_t"], theta["blu_l2_t"], theta["blu_l3_t"], theta["blu_l4_t"], theta["blu_l5_t"], theta["blu_l6_t"], theta["blu_l7_t"],  theta["blu_l7_t"], elem_window, elem_glass_blu]
	
	#update other red spect elem throughputs:
	#0 is collimator
	#1 is DichroicRG
	for i,elem in enumerate(red_elem_list):
		e = 2 + i
		for surf in llamas_red.elements[e].surfaces:
			surf.setThroughputTab(elem.fine_domain, thru(elem))

	#update other green spect elem throughputs:
	#0 is collimator
	#1 is DichroicRG
	#2 is DichroicBG
	for i,elem in enumerate(gre_elem_list):
		e = 3 + i
		for surf in llamas_green.elements[e].surfaces:
			surf.setThroughputTab(elem.fine_domain, thru(elem))
			
	#update other blue spect elem throughputs:
	#0 is collimator
	#1 is DichroicRG
	#2 is DichroicBG
	for i,elem in enumerate(blu_elem_list):
		e = 3 + i
		for surf in llamas_blue.elements[e].surfaces:
			surf.setThroughputTab(elem.fine_domain, thru(elem))
	
	#update sensors
	for llamas_channel in [llamas_red, llamas_green, llamas_blue]:
		llamas_channel.sensor.naxis1 = nx
		llamas_channel.sensor.naxis2 = ny
		llamas_channel.sensor.pixelsize = pixel_size_mm
	
	llamas_red.sensor.rn = rn_red
	llamas_red.sensor.dark = dc_red
	llamas_red.sensor.qe.setThroughputTab(elem_qe_red.fine_domain, thru(elem_qe_red))
	
	llamas_green.sensor.rn = rn_gre
	llamas_green.sensor.dark = dc_gre
	llamas_green.sensor.qe.setThroughputTab(elem_qe_gre.fine_domain, thru(elem_qe_gre))
	
	llamas_blue.sensor.rn = rn_blu
	llamas_blue.sensor.dark = dc_blu
	llamas_blue.sensor.qe.setThroughputTab(elem_qe_blu.fine_domain, thru(elem_qe_blu))
	
	#And at last, we recalculate the throughputs based on the new data
	llamas_red.throughput = llamas_red.calc_throughput(waves)
	llamas_green.throughput = llamas_green.calc_throughput(waves)
	llamas_blue.throughput = llamas_blue.calc_throughput(waves)
	
	if verbose:
		llamas_plot_throughput.plot_throughput(llamas_red, llamas_green, llamas_blue)
	#breakpoint()
	
	##############################
	###Run model
	##############################
	
	"""
	#This uses the Emission Sensitivity requirement & model
	# LSB emission 3e-19 erg/cm2/s/A/fiber SNR=5 in 1 night or less
	#simulated with an oversampled spectrum
	
	sampling = 0.05 
	wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
	#source = 3e-19 #from thesis
	#source = 5.3e-18 #dont remember
	#source = 1.71e-17 #dont remember
	source = 4.5e-18 #This is from Rob's 7/14/22 email
	spec_oversampled = wave_oversampled * 0 + source
	
	red_photons, red_noise = observe.observe_spectrum(llamas_red, tau, wave_oversampled, spec_oversampled, skyspec)
	gre_photons, gre_noise = observe.observe_spectrum(llamas_green, tau, wave_oversampled, spec_oversampled, skyspec)
	blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, tau, wave_oversampled, spec_oversampled, skyspec)
	spectral_npix = 4
	signal = np.array(spectral_npix * (blu_photons + gre_photons + red_photons))
	noise = np.array(np.sqrt(spectral_npix) * (blu_noise + gre_noise + red_noise))
	snr = [s/n for s,n in zip(signal, noise)]
	mean_snr = np.mean(snr)
	"""
	
	#This uses the Continuum Sensitivity requirement & model
	# R = 25 Lyman Break Galaxy SNR=3 in 5 hours
	#simulated with a Shapley spectrum
	input_wv = x["shapley_wv"]
	input_spec = x["shapley_spec"]
	
	fnu_shapley = 3.4e-30
	Rmag = 25.0
	c = 3.0e17   # nm/sec
	fnu = 10**(-0.4*(Rmag+48.6))
	input_spec *= fnu / fnu_shapley
	input_flam = input_spec * c / (input_wv)**2 
	
	red_photons, red_noise = observe.observe_spectrum(llamas_red, tau, input_wv, input_flam, skyspec)
	gre_photons, gre_noise = observe.observe_spectrum(llamas_green, tau, input_wv, input_flam, skyspec)
	blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, tau, input_wv, input_flam, skyspec)
	signal = np.array(blu_photons + gre_photons + red_photons)
	noise = np.array(blu_noise + gre_noise + red_noise)
	snr = [s/n for s,n in zip(signal, noise)]
	mean_snr = np.mean(snr)
	
	"""
	for i,s in enumerate(signal):
		if s == 0:
			print(waves[i], "signal zero")
		if math.isnan(s):
			print(waves[i], "signal nan")

	for i,s in enumerate(noise):
		if s == 0:
			print(waves[i], "signal zero")
		if math.isnan(s):
			print(waves[i], "noise nan")
	"""
			
	return mean_snr
	

