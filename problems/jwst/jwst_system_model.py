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
		matlab.double(x["caxi"]),
		matlab.double(x["crad"]),
		matlab.double(x["ctor"]),
		matlab.double(x["haxi"]),
		matlab.double(x["hrad"]),
		matlab.double(x["htor"]),
		matlab.double(x["zeta1"]),
		matlab.double(x["a"]),
		matlab.double(x["wheel_locs"]),
		matlab.double(x["n"]),
		matlab.double(x["h"]),
		matlab.double(x["C"]),
		matlab.double(x["Nsurf"]),
		matlab.double(x["D"]),
		matlab.double(x["BP"]),
		matlab.double(x["PH"]),
		matlab.double(x["R0"]),
		matlab.double(x["mass"]),
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
		1 if verbose else 0 #diagnostics
	)
	
	return jitter
	
	
if __name__ == '__main__':  
	#This ought to be done at the problem level	
	0