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

#system model startup -- is it ok to have one universal instance?
import matlab.engine

#setup matlab engine
eng = matlab.engine.start_matlab()
eng.cd(r'nexus_lyapunov', nargout=0)

#Wrapper function for model call
#def sensitivity_hlva(x, d, col, config_table, spectrograph_models):
def nexus_lyapunov_system_model(theta, x, verbose=True):
	#call matlab function
	jitter = eng.nexus_lyapunov(
		matlab.double(x["Ru"]),
		matlab.double(theta["Us"]),
		matlab.double(theta["Ud"]),
		matlab.double(x["fc"]),
		matlab.double(theta["Qc"]),
		matlab.double(x["Tst"]),
		matlab.double(x["Srg"]),
		matlab.double(x["Sst"]),
		matlab.double(x["Tgs"]),
		matlab.double(x["m_SM"]),
		matlab.double(x["m_SMhub"]),
		matlab.double(x["I_SMhubt"]),
		matlab.double(x["I_SMhuba"]),
		matlab.double(x["m_RW"]),
		matlab.double(theta["K_yPM"]),
		matlab.double(x["I_xRWA"]),
		matlab.double(x["I_yRWA"]),
		matlab.double(x["I_RWt"]),
		matlab.double(x["I_RWa"]),
		matlab.double(x["m_ISO"]),
		matlab.double(x["I_ISOa"]),
		matlab.double(x["I_ISOt"]),
		matlab.double(theta["K_yISO"]),
		matlab.double(theta["K_xISO"]),
		matlab.double(x["m_RWAchx"]),
		matlab.double(x["I_bus"]),
		matlab.double(x["m_bus"]),
		matlab.double(x["m_prop"]),
		matlab.double(x["I_propt"]),
		matlab.double(x["I_propa"]),
		matlab.double(x["m_instr"]),
		matlab.double(x["I_i1"]),
		matlab.double(x["I_i2"]),
		matlab.double(x["I_i3"]),
		matlab.double(x["A_sptop"]),
		matlab.double(x["D_sp"]),
		matlab.double(x["t_sp"]),
		matlab.double(x["I_ss"]),
		matlab.double(theta["K_rad1"]),
		matlab.double(theta["K_rad2"]),
		matlab.double(theta["K_rISO"]),
		matlab.double(theta["K_act1"]),
		matlab.double(theta["K_act2"]),
		matlab.double(x["I_iso"]),
		matlab.double(theta["K_zpet"]),
		matlab.double(x["lambda_"]),
		matlab.double(x["Ro"]),
		matlab.double(x["QE"]),
		matlab.double(x["Mgs"]),
		matlab.double(x["fca"]),
		matlab.double(x["Kc"]),
		matlab.double(x["Kcf"]),
		matlab.double(x["nray"]),
		matlab.double([2.95000000000000e-08,1.41600000000000e-08,7.14000000000000e-09, .36000000000000e-09,2.49100000000000e-08,1.08700000000000e-08, 3.01000000000000e-08,1.36200000000000e-08,7.05000000000000e-09]),
		matlab.double([7.851999999999999e-08,2.984000000000000e-08,7.770000000000000e-09, 1.036000000000000e-08, 1.280000000000000e-08, 8.490000000000000e-09,9.850000000000000e-09, 1.849000000000000e-08, 6.180000000000000e-09, 4.020000000000000e-09]),
		matlab.double([3.239000000000000e-08,7.439999999999999e-09,3.090000000000000e-09, 4.719999999999999e-09,1.210000000000000e-09,4.120000000000000e-09, 9.239999999999999e-09,1.027000000000000e-08,6.300000000000000e-09, 8.510000000000000e-09]),
		matlab.double([1,2,2.900000000000000,3.880000000000000,4.290000000000000,4.850000000000000,5.400000000000000,5.810000000000000,6.160000000000000]),
		matlab.double([1,2,3,4,4.430000000000000,5,5.380000000000000,5.600000000000000,6, 6.390000000000000]),
		matlab.double([1,2,2.900000000000000,3.880000000000000,4,4.430000000000000, 5.200000000000000,5.400000000000000,5.600000000000000,5.810000000000000]),
		matlab.double(x["zeta1"]),
		matlab.double(x["a"]),
		matlab.double([	  
			[79.00000000000000,-0.47413000000000,-0.81375000000000, 2.04467000000000],
			[80.00000000000000,-0.47413000000000,-0.96331000000000, 1.72393000000000],
			[81.00000000000000,-0.47413000000000,-1.28405000000000, 1.87349000000000],
			[82.00000000000000,-0.47413000000000,-1.13448000000000, 2.19423000000000],
			[83.00000000000000,-0.56728000000000,-1.04890000000000, 1.95908000000000]
		]),
		matlab.double(x["n"]),
		matlab.double([1.0,2.0,3.0,4.0,5.0,6.0]),
		matlab.double([
			[42, 0.95,4.1, 2.75,0.9,1.2],
			[0.2,0.09,0.25,1.0, 5.0,0.4]
		]),
		matlab.double(x["Nsurf"]),
		matlab.double(x["D"]),
		matlab.double(x["BP"]),
		matlab.double(x["PH"]),
		matlab.double(x["R0"]),
		matlab.double([
			[752.094612601836,1.03216046820620e-15,2.23432383705813e-15,3.23726195090639e-17,1.47160061914065e-13,4.28546087505310e-14],
			[5.79397640976254e-16,752.094612601836,-4.27435864480685e-14,-1.04638520070921e-13,2.43945523111095e-15,-6.75674793892966e-15],
			[2.50494069931051e-15,-4.03219124756049e-14,752.094612601836,-1.60062087509864e-13,4.98039109952941e-15,-1.69109226052856e-15],
			[-7.47073394395928e-16,-1.47049039611602e-13,-1.31672450720544e-13,1720.96055133990,-10.5270557086614,13.3850421894362],
			[1.47160061914065e-13,2.62130691004106e-15,4.98039109952941e-15,-10.5270557086614,1352.56559956577,228.801871490305],
			[4.28546087505310e-14,-6.75674793892966e-15,-2.65496910573658e-15,13.3850421894362,228.801871490305,1508.50841532254]
		]),
		matlab.double(x["FgsNom"]),
		matlab.double(theta["K_pm1"]),
		matlab.double(theta["K_pm3"]),
		matlab.double(theta["K_pm4"]),
		matlab.double(theta["K_pm5"]),
		matlab.double(theta["K_pm6"]),
		matlab.double(theta["K_act_pm2"]),
		matlab.double(theta["K_act_pm3"]),
		matlab.double(theta["K_act_pm4"]),
		matlab.double(theta["K_act_pm5"]),
		matlab.double(theta["K_act_pm6"]),
		matlab.double(theta["K_xpet"]),
		matlab.double(theta["c_RWA"]),
		matlab.double(theta["c_RWAI"]),
		matlab.double(theta["c_SM_act"]),
		matlab.double(theta["c_PM"]),
		matlab.double(theta["c_PM_act"]),
		matlab.double(theta["c_petal"]),
		matlab.double(theta["zeta_sunshield"]),
		matlab.double(theta["zeta_isolator"]),
		matlab.double(theta["zeta_solarpanel"]),
		matlab.double(theta["CCA_factor"]),
		1 if verbose else 0 #diagnostics
	)
	
	return jitter
	
	
if __name__ == '__main__':  
	#This ought to be done at the problem level	
	0