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

from problems.llamas_snr import *

"""
Full matrix experiment model
"""
def fp_likelihood_fn(theta, d, x, err=True):
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