#trying to implement the ideas in Lieberman & Willcos 2014

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

import sys.path.append('../..')
#focal plane
from problems.fp_verification.fp_problem import *

d_historical = [
				20,   #t_gain
				30,   #I_gain
				1,	#n_meas_rn
				8,	#d_num
				9600, #d_max
				2	 #d_pow   #approx
			   ]

def expt_process(theta):
	yd = fp.eta(theta, d_historical)
	return yd

def pred_process(theta):
	yp = fp.H(theta)
	return yp
	
def theta_prior_sample(num):
	return fp.prior_rvs(num)
	
def theta_prior_pdf(theta):
	return fp.prior_pdf_unnorm(theta)
	
