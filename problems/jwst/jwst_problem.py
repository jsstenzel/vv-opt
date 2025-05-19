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
#from problems.jwst.jwst_cost_model import *
import matlab.engine

def unif_margin(nominal):
	return [nominal - 0.2*nominal, nominal + 0.2*nominal]

def construct_jwst_jitter_problem(verbose_probdef=False):
	print("Constructing jwst problem ...",flush=True)

	theta_defs = [
					["Us", ["uniform", unif_margin(1.8)], "continuous"],
					["Ud", ["uniform", unif_margin(60.0)], "continuous"],
					["Qc", ["uniform", unif_margin(0.005)], "continuous"],
					["I_SMhubt", ["uniform", unif_margin(0.25200E+00)], "continuous"],
					["I_SMhuba", ["uniform", unif_margin(0.45900E+00)], "continuous"],
					["K_yPM", ["uniform", unif_margin(0.77400E+06)], "continuous"],
					["I_xRWA", ["uniform", unif_margin(0.40187E+00)], "continuous"],
					["I_yRWA", ["uniform", unif_margin(0.22445E+00)], "continuous"],
					["I_RWt", ["uniform", unif_margin(0.83595E-01)], "continuous"],
					["I_RWa", ["uniform", unif_margin(0.14339E+00)], "continuous"],
					["I_ISOa", ["uniform", unif_margin(0.11720E-03)], "continuous"],
					["I_ISOt", ["uniform", unif_margin(0.39300E-01)], "continuous"],
					["K_yISO", ["uniform", unif_margin(0.14600E+04)], "continuous"],
					["K_xISO", ["uniform", unif_margin(0.14000E+12)], "continuous"],
					["I_bus", ["uniform", unif_margin(0.85080E+02)], "continuous"],
					["I_propt", ["uniform", unif_margin(0.51100E+01)], "continuous"],
					["I_propa", ["uniform", unif_margin(0.74000E+00)], "continuous"],
					["I_i1", ["uniform", unif_margin(0.49200E+01)], "continuous"],
					["I_i2", ["uniform", unif_margin(0.75420E+01)], "continuous"],
					["I_i3", ["uniform", unif_margin(0.41280E+01)], "continuous"],
					["A_sptop", ["uniform", unif_margin(0.14040E-02)], "continuous"],
					["D_sp", ["uniform", unif_margin(0.060)], "continuous"],
					["t_sp", ["uniform", unif_margin(0.003)], "continuous"],
					["I_ss", ["uniform", unif_margin(0.78350E-08)], "continuous"],
					["K_rad1", ["uniform", unif_margin(0.50000E+6)], "continuous"],
					["K_rad2", ["uniform", unif_margin(0.30000E+6)], "continuous"],
					["K_rISO", ["uniform", unif_margin(3000)], "continuous"],
					["K_act1", ["uniform", unif_margin(0.20000E+11)], "continuous"],
					["K_act2", ["uniform", unif_margin(0.14000E+12)], "continuous"],
					["I_iso", ["uniform", unif_margin(1.00000E-5)], "continuous"],
					["K_zpet", ["uniform", unif_margin(0.9000E+08)], "continuous"],
					["K_pm1", ["uniform", unif_margin(0.10000E+07)], "continuous"],
					["K_pm3", ["uniform", unif_margin(0.58400E+06)], "continuous"],
					["K_pm4", ["uniform", unif_margin(0.59820E+02)], "continuous"],
					["K_pm5", ["uniform", unif_margin(0.49000E+02)], "continuous"],
					["K_pm6", ["uniform", unif_margin(0.33250E+02)], "continuous"],
					["K_act_pm2", ["uniform", unif_margin(0.29100E+07)], "continuous"],
					["K_act_pm3", ["uniform", unif_margin(0.10000E+07)], "continuous"],
					["K_act_pm4", ["uniform", unif_margin(0.33250E+02)], "continuous"],
					["K_act_pm5", ["uniform", unif_margin(0.49000E+02)], "continuous"],
					["K_act_pm6", ["uniform", unif_margin(0.12012E+03)], "continuous"],
					["K_xpet", ["uniform", unif_margin(1e16)], "continuous"],
					["c_RWA", ["uniform", unif_margin(0.01*math.sqrt(0.25000E+01)*math.sqrt(0.14600E+04))], "continuous"],
					["c_RWAI", ["uniform", unif_margin(0.01*math.sqrt(0.15E+01)*math.sqrt(0.14600E+04))], "continuous"],
					["c_SM_act", ["uniform", unif_margin(0.01*math.sqrt(2.49)*math.sqrt(0.30000E+6))], "continuous"],
					["c_PM", ["uniform", unif_margin(0.01*math.sqrt(0.18860E+02)*math.sqrt(0.77400E+06))], "continuous"],
					["c_PM_act", ["uniform", unif_margin(0.01*math.sqrt(0.18860E+02)*math.sqrt(0.30000E+6))], "continuous"],
					["c_petal", ["uniform", unif_margin(0.01*math.sqrt(0.18860E+02)*math.sqrt(0.9000E+08))], "continuous"],
					["zeta_sunshield", ["uniform", unif_margin(0.005*5)], "continuous"],
					["zeta_isolator", ["uniform", unif_margin(0.005*20)], "continuous"],
					["zeta_solarpanel", ["uniform", unif_margin(0.005*20)], "continuous"],
				]

	y_defs = [	
					"y_gain_red", 
				]

	d_defs = [
					["t_gain", ['uniform', [.1, 600]], "continuous"], #gain
				]

	x_defs = [
				#put eng in here, for accessibility
				#["eng", [], "discrete", eng],
				#x dict fixed
				["m_SM", [], "discrete", 2.49],
				["m_SMhub", [], "discrete", 4.4],
				["m_RW", [], "discrete", 0.10600E+02],
				["m_bus", [], "discrete", 0.33200E+03],
				["m_prop", [], "discrete", 0.37000E+02],
				["m_ISO", [], "discrete", 0.15E+01],
				["m_RWAchx", [], "discrete", 0.25000E+01],
				["m_instr", [], "discrete", 0.50000E+02],
				["wheel_locs", [], "discrete", [	  
					[79.00000000000000,-0.47413000000000,-0.81375000000000, 2.04467000000000],
					[80.00000000000000,-0.47413000000000,-0.96331000000000, 1.72393000000000],
					[81.00000000000000,-0.47413000000000,-1.28405000000000, 1.87349000000000],
					[82.00000000000000,-0.47413000000000,-1.13448000000000, 2.19423000000000],
					[83.00000000000000,-0.56728000000000,-1.04890000000000, 1.95908000000000]
				]],
				["nray", [], "discrete",  134],
				["caxi", [], "discrete",  [2.95000000000000e-08,1.41600000000000e-08,7.14000000000000e-09, .36000000000000e-09,2.49100000000000e-08,1.08700000000000e-08, 3.01000000000000e-08,1.36200000000000e-08,7.05000000000000e-09]],
				["crad", [], "discrete",  [7.851999999999999e-08,2.984000000000000e-08,7.770000000000000e-09, 1.036000000000000e-08, 1.280000000000000e-08, 8.490000000000000e-09,9.850000000000000e-09, 1.849000000000000e-08, 6.180000000000000e-09, 4.020000000000000e-09]],
				["ctor", [], "discrete",  [3.239000000000000e-08,7.439999999999999e-09,3.090000000000000e-09, 4.719999999999999e-09,1.210000000000000e-09,4.120000000000000e-09, 9.239999999999999e-09,1.027000000000000e-08,6.300000000000000e-09, 8.510000000000000e-09]],
				["haxi", [], "discrete",  [1,2,2.900000000000000,3.880000000000000,4.290000000000000,4.850000000000000,5.400000000000000,5.810000000000000,6.160000000000000]],
				["hrad", [], "discrete",  [1,2,3,4,4.430000000000000,5,5.380000000000000,5.600000000000000,6, 6.390000000000000]],
				["htor", [], "discrete",  [1,2,2.900000000000000,3.880000000000000,4,4.430000000000000, 5.200000000000000,5.400000000000000,5.600000000000000,5.810000000000000]],
				["zeta1", [], "discrete", 0.2],
				["a", [], "discrete", 0.2],
				["n", [], "discrete", 3],
				["h", [], "discrete", [1.0,2.0,3.0,4.0,5.0,6.0]],
				["C", [], "discrete", [
					[42,0.95,4.1,2.75,0.9,1.2],
					[0.2,0.09,0.25,1.0,5.0,0.4]
				]],
				["Nsurf", [], "discrete", 10],
				["D", [], "discrete", 2.8],
				["BP", [], "discrete", 0.4],
				["PH", [], "discrete", 1.0e10],
				["R0", [], "discrete", 4*60],
				["mass", [], "discrete", [
					[752.094612601836,1.03216046820620e-15,2.23432383705813e-15,3.23726195090639e-17,1.47160061914065e-13,4.28546087505310e-14],
					[5.79397640976254e-16,752.094612601836,-4.27435864480685e-14,-1.04638520070921e-13,2.43945523111095e-15,-6.75674793892966e-15],
					[2.50494069931051e-15,-4.03219124756049e-14,752.094612601836,-1.60062087509864e-13,4.98039109952941e-15,-1.69109226052856e-15],
					[-7.47073394395928e-16,-1.47049039611602e-13,-1.31672450720544e-13,1720.96055133990,-10.5270557086614,13.3850421894362],
					[1.47160061914065e-13,2.62130691004106e-15,4.98039109952941e-15,-10.5270557086614,1352.56559956577,228.801871490305],
					[4.28546087505310e-14,-6.75674793892966e-15,-2.65496910573658e-15,13.3850421894362,228.801871490305,1508.50841532254]
				]],
				["FgsNom", [], "discrete", 30.0],
				
				#x dict random
				["Ru", ["uniform", unif_margin(3000)], "continuous", 3000.0],
				["fc", ["uniform", unif_margin(30)], "continuous", 30.0],
				["Tst", ["uniform", unif_margin(20)], "continuous", 20.0],
				["Srg", ["uniform", unif_margin(3e-14)], "continuous", 3e-14],
				["Sst", ["uniform", unif_margin(2)], "continuous", 2],
				["Tgs", ["uniform", unif_margin(0.04)], "continuous", 0.04],
				["lambda_", ["uniform", unif_margin(1e-6)], "continuous", 1e-6],
				["Ro", ["uniform", unif_margin(0.98)], "continuous", 0.98],
				["QE", ["uniform", unif_margin(0.8)], "continuous", 0.8],
				["Mgs", ["uniform", unif_margin(15)], "continuous", 15],
				["fca", ["uniform", unif_margin(0.01)], "continuous", 0.01],
				["Kc", ["uniform", [0,1]], "continuous", 0.0],
				["Kcf", ["uniform", unif_margin(2000)], "continuous", 2000.0],
			]

	#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
	eta = jwst_eta
	H = nexus_lyapunov_system_model
	Gamma = None#jwst_cost
	jwst_jitter = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs, x_defs)
	print("jwst_jitter problem constructed.",flush=True)
	return jwst_jitter
	
if __name__ == '__main__':  
	jwst_jitter = construct_jwst_jitter_problem(verbose_probdef=True)
	print(jwst_jitter)
	
"""
theta nominal:
	Us=1.8
	Ud=60.0
	Qc=0.005
	I_SMhubt=0.25200E+00
	I_SMhuba=0.45900E+00
	K_yPM=0.77400E+06
	I_xRWA=0.40187E+00
	I_yRWA=0.22445E+00
	I_RWt=0.83595E-01
	I_RWa=0.14339E+00
	I_ISOa=0.11720E-03
	I_ISOt=0.39300E-01
	K_yISO=0.14600E+04
	K_xISO=0.14000E+12
	I_bus=0.85080E+02
	I_propt=0.51100E+01
	I_propa=0.74000E+00
	I_i1=0.49200E+01
	I_i2=0.75420E+01
	I_i3=0.41280E+01
	A_sptop=0.14040E-02
	D_sp=0.060
	t_sp=0.003
	I_ss=0.78350E-08
	K_rad1=0.50000E+6
	K_rad2=0.30000E+6
	K_rISO=3000
	K_act1=0.20000E+11
	K_act2=0.14000E+12
	I_iso=1.00000E-5
	K_zpet=0.9000E+08
	K_pm1=0.10000E+07
	K_pm3=0.58400E+06
	K_pm4=0.59820E+02
	K_pm5=0.49000E+02
	K_pm6=0.33250E+02
	K_act_pm2=0.29100E+07
	K_act_pm3=0.10000E+07
	K_act_pm4=0.33250E+02
	K_act_pm5=0.49000E+02
	K_act_pm6=0.12012E+03
	K_xpet=1e16
	c_RWA=0.01*math.sqrt(0.25000E+01)*math.sqrt(0.14600E+04)
	c_RWAI=0.01*math.sqrt(0.15E+01)*math.sqrt(0.14600E+04)
	c_SM_act=0.01*math.sqrt(2.49)*math.sqrt(0.30000E+6)
	c_PM=0.01*math.sqrt(0.18860E+02)*math.sqrt(0.77400E+06)
	c_PM_act=0.01*math.sqrt(0.18860E+02)*math.sqrt(0.30000E+6)
	c_petal=0.01*math.sqrt(0.18860E+02)*math.sqrt(0.9000E+08)
	zeta_sunshield=0.005*5
	zeta_isolator=0.005*20
	zeta_solarpanel=0.005*20
"""