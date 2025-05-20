import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv
import dill

sys.path.append('../..')
#jwst
from problems.jwst.jwst_problem import *
#analysis
from obed.obed_multivar import *
from obed.obed_gbi import *
#from obed.pdf_estimation import *
from inference.bn_modeling import *
from uq.uncertainty_propagation import *
from uq.sensitivity_analysis import *
from uq.saltelli_gsa import *
from uq.gsa_convergence import *
from opt.ngsa import *

from scipy.optimize import minimize

problem = construct_jwst_jitter_problem()

def damping_to_jitter(problem, c_RWA, c_RWAI, c_SM_act, c_PM, c_PM_act, c_petal):
	theta_augment = problem.theta_nominal
	theta_augment[42] = c_RWA
	theta_augment[43] = c_RWAI
	theta_augment[44] = c_SM_act
	theta_augment[45] = c_PM
	theta_augment[46] = c_PM_act
	theta_augment[47] = c_petal
	
	QoI_nominal = problem.H(theta_augment, verbose=True)
	print("Given the nominal theta:", theta_nominal)
	#print("Nominal y:", y_nominal)
	print("Nominal QoI:", QoI_nominal)
	
def damping_scale_to_jitter(problem, S, verbose=True):
	theta_augment = problem.theta_nominal
	theta_augment[42] *= S
	theta_augment[43] *= S
	theta_augment[44] *= S
	theta_augment[45] *= S
	theta_augment[46] *= S
	theta_augment[47] *= S
	
	QoI = problem.H(theta_augment, verbose=verbose)
	return QoI
	
def fn_to_minimize(S):
	return abs(15.5134 - damping_scale_to_jitter(problem, S, verbose=False))

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Function to run for this vvopt analysis')
	parser.add_argument('-n', type=float, default=0, help='Number of iterations to give to the function')
	parser.add_argument('--filename', metavar='string', default="SA_QoI", help='Base name to five to SA_QoI_sample')
	args = parser.parse_args()

	if args.run == "single_component_damping":
		damping = float(args.n)
		jitter = damping_to_jitter(problem, c_RWA=damping, c_RWAI=damping, c_SM_act=damping, c_PM=damping, c_PM_act=damping, c_petal=damping)
		print("Jitter:", jitter)
	
	elif args.run == "scale_damping":
		scale = float(args.n)
		jitter = damping_scale_to_jitter(problem, scale)
		print("Jitter:", jitter)
		
	elif args.run == "train_scale_damping":
		x0 = 1e-6
		bounds = [(0,None)]
		res = minimize(fn_to_minimize, x0, method='Nelder-Mead', tol=1e-6, bounds=bounds)
		print(res.x)
		print(damping_scale_to_jitter(problem, res.x[0]))
		print(fn_to_minimize(res.x[0]))
			
	#right scale is between 0.000001 and 0.0000001
	#between 4e-7 (15.07) and 4.5e-7 (15.93)
	#4.23245e-7 is good enough
	
	else:
		print("I dont recognize the command",args.run)
