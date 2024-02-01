#import argparse
import os
import sys
import shutil
import csv
import fileinput
sys.path.insert(0, "..")

import spectrograph as spec
import matplotlib.pyplot as plt
import observe
from astropy.io import fits
import scipy.signal as ss
import spectrograph as spec
import observe

import numpy as np
import itertools
import multiprocessing as mp
import math
from SALib.sample import saltelli, fast_sampler
from SALib.analyze import sobol, fast
from copy import deepcopy
import scipy.optimize as optimization

os.environ['COATINGS_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/"
os.environ['TEST_DATA_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/TEST DATA/"

num_workers = mp.cpu_count() 

################################
#Utility functions
################################

def mkdir(path):
	if os.path.exists(path):
		shutil.rmtree(path)
	os.makedirs(path)
	
def soft_mkdir(path):
	if os.path.exists(path):
		return
	os.makedirs(path)
	
def update_tab(thru_tab, bias):
	return_thru = []
	
	for thru in thru_tab:
		new_thru = thru + bias
		if new_thru > 1:
			new_thru = 1
		if new_thru < 0:
			new_thru = 0
		return_thru.append(new_thru)
	return_thru = np.array(return_thru)
	return_refl = 1 - return_thru
	return [return_thru, return_refl]
	
def parabola(x,a,b,c):
    return a*x*x + b*x + c	

def make_vph_tab(points,errors,color):
	thrus = [float(p)+float(e) for p,e in zip(points,errors)]
	waves = []
	if color == "r":
		waves = [700.0,825.0,925.0]
	if color == "g":
		waves = [475.0,550.0,675.0]
	if color == "b":
		waves = [375.0,425.0,475.0]
	fit = optimization.curve_fit(parabola,waves,thrus,[-.000005,.005,.5])
	out_waves = np.arange(350,975,.5)
	out_trans = [min(1,max(0,parabola(wave,fit[0][0],fit[0][1],fit[0][2]))) for wave in out_waves]
	out_reflect = [1-t for t in out_trans]
	
	return [out_trans, out_reflect, out_waves]


################################
#Utility functions
################################

#Wrapper function for model call
def sensitivity_hlva(x, d, col, config_table, spectrograph_models, skyspec):
	spect_models = deepcopy(spectrograph_models)
    
    #Expected format for x: 
    x_names = ["resolution","f_col","f_cam","fratio_col","fratio_cam","blade_obscure","d_beam","fiber_theta","fiber_core_d","fiber_frd","microlens_thru","fiber_ar_thru","fiber_internal_thru","collimator_thru","sl_thru","bg_thru","red_vph","green_vph","blue_vph","red_prism","green_prism","blue_prism","l1_red_thru","l2_red_thru","l3_red_thru","l4_red_thru","l5_red_thru","l6_red_thru","l7_red_thru","red_window_thru","l1_green_thru","l2_green_thru","l3_green_thru","l4_green_thru","l5_green_thru","l6_green_thru","l7_green_thru","green_window_thru","l1_blue_thru","l2_blue_thru","l3_blue_thru","l4_blue_thru","l5_blue_thru","l6_blue_thru","l7_blue_thru","blue_window_thru","red_internal_thru","green_internal_thru","blue_internal_thru","red_rn","red_dc","green_rn","green_dc","blue_rn","blue_dc","skyfile"]
]
    make a map

    #Expected format for d:
    d_names = 
    make a map
    
    ###Take in copied x values
    x_u = deepcopy(x)
    
    ###Simulate experiments w/ d to find uncertainties
    ###Apply uncertainties to x
    exp_
    
    ###Load in the new values to the spect_models
    """
	#update frd
	llamas_red.fiber.frd_loss += frd_error
	llamas_green.fiber.frd_loss += frd_error
	llamas_blue.fiber.frd_loss += frd_error
	#update cameras
	llamas_red.sensor.rn += red_rn_error
	llamas_red.sensor.dark += red_dc_error
	llamas_green.sensor.rn += gre_rn_error
	llamas_green.sensor.dark += gre_dc_error
	llamas_blue.sensor.rn += blu_rn_error
	llamas_blue.sensor.dark += blu_dc_error
	#update sl
	llamas_red.elements[1].surfaces[0].setThroughputTab(update_tab(llamas_red.elements[1].surfaces[0].transmission_tab,sl_error)) 
	llamas_green.elements[1].surfaces[0].setThroughputTab(update_tab(llamas_green.elements[1].surfaces[0].transmission_tab,-sl_error)) 
	llamas_blue.elements[1].surfaces[0].setThroughputTab([llamas_green.elements[1].surfaces[0].transmission_tab,llamas_green.elements[1].surfaces[0].reflectance_tab]) 
	#update bg
	llamas_green.elements[2].surfaces[0].setThroughputTab(update_tab(llamas_green.elements[2].surfaces[0].transmission_tab,bg_error)) 
	llamas_blue.elements[2].surfaces[0].setThroughputTab(update_tab(llamas_blue.elements[2].surfaces[0].transmission_tab,-bg_error))
	#update vph
	red_vph_points = [float(row[col['red_vph_pt1']]),float(row[col['red_vph_pt2']]),float(row[col['red_vph_pt3']])]  
	green_vph_points = [float(row[col['green_vph_pt1']]),(row[col['green_vph_pt2']]),float(row[col['green_vph_pt3']])]  
	blue_vph_points = [float(row[col['blue_vph_pt1']]),float(row[col['blue_vph_pt2']]),float(row[col['blue_vph_pt3']])]  
	red_vph_tabs = make_vph_tab(red_vph_points,[red_vph_err1,red_vph_err2,red_vph_err3],'r')
	green_vph_tabs = make_vph_tab(green_vph_points,[green_vph_err1,green_vph_err2,green_vph_err3],'r')
	blue_vph_tabs = make_vph_tab(blue_vph_points,[blue_vph_err1,blue_vph_err2,blue_vph_err3],'r')
	llamas_red.grating.blaze.setThroughputTab(red_vph_tabs)
	llamas_green.grating.blaze.setThroughputTab(green_vph_tabs)
	llamas_blue.grating.blaze.setThroughputTab(blue_vph_tabs)
	#update lens throughputs (spect 12345678, camera rgb, lens L1 L2 L3 L4-L5 L6-L7 W, two surfaces)
	set up iterator through lens_rgb_...
	for lens in [4,5,6,7,8,9,10,11]:
		for surf in [0,1]:
			llamas_red.elements[lens-1].surfaces[surf].setThroughputTab(update_tab(llamas_red.elements[lens-1].surfaces[surf].transmission_tab,lens_thru_err)) 
			llamas_green.elements[lens].surfaces[surf].setThroughputTab(update_tab(llamas_red.elements[lens].surfaces[surf].transmission_tab,lens_thru_err)) 
			llamas_blue.elements[lens].surfaces[surf].setThroughputTab(update_tab(llamas_red.elements[lens].surfaces[surf].transmission_tab,lens_thru_err)) 
	#update fiber core
	llamas_red.fiber.dFib += fiber_core_d_err
	llamas_green.fiber.dFib += fiber_core_d_err
	llamas_blue.fiber.dFib += fiber_core_d_err
    """

	llamas_red.throughput = llamas_red.calc_throughput(waves)
	llamas_green.throughput = llamas_green.calc_throughput(waves)
	llamas_blue.throughput = llamas_blue.calc_throughput(waves)
    
    ###Run model
	run_snrs = []
	
	texp=600
	sampling = 0.05 
	wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
	spec_oversampled = wave_oversampled * 0 + 1.71e-17 #5.3e-18
	
	red_photons, red_noise = observe.observe_spectrum(llamas_red, texp, wave_oversampled, spec_oversampled, skyspec)
	gre_photons, gre_noise = observe.observe_spectrum(llamas_green, texp, wave_oversampled, spec_oversampled, skyspec)
	blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, texp, wave_oversampled, spec_oversampled, skyspec)
	spectral_npix = 4
	signal = np.array(spectral_npix * (blu_photons + gre_photons + red_photons))
	noise = np.array(np.sqrt(spectral_npix) * (blu_noise + gre_noise + red_noise))
	snr = [signal[i]/noise[i] for i in range(0, len(signal))]
	run_snrs.append(np.median(snr))
			
	return(np.mean(run_snrs))
	

