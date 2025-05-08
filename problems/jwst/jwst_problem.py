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
from problems.jwst.jwst_exp_models import *
from problems.jwst.jwst_system_model import *
from problems.jwst.jwst_cost_model import *


def construct_llamas_snr_problem(verbose_probdef=False):
	print("Constructing jwst priors ...",flush=True)

	#define priors
	prior_gain_SN1 = ["gamma_mv",  [0.999,0.2**2]] #mean, variance

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
					["ny", ["nonrandom", [2048]], "discrete", 2048],
					["pixel_size_mm", [], "continuous", 13.5],
				]

	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	eta = jwst_eta
	H = jwst
	Gamma = snr_cost
	llamas_snr = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs, x_defs)
	print("llamas_snr_full problem constructed.",flush=True)
	return llamas_snr
	
if __name__ == '__main__':  
	llamas_snr = construct_llamas_snr_problem(verbose_probdef=True)
	print(llamas_snr)