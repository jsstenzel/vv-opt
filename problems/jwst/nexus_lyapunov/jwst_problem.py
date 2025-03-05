import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv

sys.path.append('../..')
from problems.problem_definition import *
from approx.gaussian_process import *
from approx.learn_gp import *
from approx.regression_models import *
#llamas
from problems.llamas_snr_full.snr_exp_models import *
from problems.llamas_snr_full.snr_system_model import *
from problems.llamas_snr_full.snr_cost_model import *

	

def construct_jwst_problem(verbose_probdef=False):
	print("Constructing llamas_snr_full priors ...",flush=True)

	#set path variables
	_dir = "./llamas-etc/COATINGS/"
	
	#Directory to load prior def files from
	sloc = './priors/'

	_wave_min = 350.0
	_wave_max = 975.0
	_wave_redgreen = 690.0
	_wave_greenblue = 480.0
	_bandpass = _wave_max - _wave_min
	#my physical argument here is that we should consider measurements at these distances apart to be essentially independent
	#and this is based on the approximate spectral resolution of the instrument
	#R = lambda / dlambda = 2200
	#so for lambda=350, dlambda=0.159
	#and for lambda=975, dlambda=0.443
	_lengthscale = 50.0

	###priors

	#these priors are based on requirements that were met, see Camera Qual Report
	theta_defs = [                             #mean, variance
						["gain_red", prior_gain_SN1, "continuous"],
					]


	y_defs = [	
					"y_gain_red", 
				]

	d_defs = [
					["t_gain", ['uniform', [.1, 600]], "continuous"], #gain
				]
		

	x_defs = [
					#focal plane
					["nx", ["nonrandom", [2048]], "discrete", 2048],
				]


	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	eta = snr_likelihood_fn
	H = sensitivity_hlva
	Gamma = snr_cost
	jwst = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs, x_defs)
	print("jwst problem constructed.",flush=True)
	return jwst
	
if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	args = parser.parse_args()
	
	jwst = construct_jwst_problem(verbose_probdef=True)
	print(jwst)