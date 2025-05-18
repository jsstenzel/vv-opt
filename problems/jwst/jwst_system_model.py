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



#Wrapper function for model call
#def sensitivity_hlva(x, d, col, config_table, spectrograph_models):
def nexus_lyapunov_system_model(theta, x, verbose=True):
	#call system model with sample theta
	#x dict fixed
	m_SM=2.49
	m_SMhub=4.4
	m_RW=0.10600E+02
	m_bus=0.33200E+03
	m_prop=0.37000E+02
	m_ISO=0.15E+01
	m_RWAchx=0.25000E+01
	m_instr=0.50000E+02
	wheel_locs=[	  
		[79.00000000000000,-0.47413000000000,-0.81375000000000, 2.04467000000000],
		[80.00000000000000,-0.47413000000000,-0.96331000000000, 1.72393000000000],
		[81.00000000000000,-0.47413000000000,-1.28405000000000, 1.87349000000000],
		[82.00000000000000,-0.47413000000000,-1.13448000000000, 2.19423000000000],
		[83.00000000000000,-0.56728000000000,-1.04890000000000, 1.95908000000000]
	]
	nray= 134
	caxi= [2.95000000000000e-08,1.41600000000000e-08,7.14000000000000e-09, .36000000000000e-09,2.49100000000000e-08,1.08700000000000e-08, 3.01000000000000e-08,1.36200000000000e-08,7.05000000000000e-09]
	crad= [7.851999999999999e-08,2.984000000000000e-08,7.770000000000000e-09, 1.036000000000000e-08, 1.280000000000000e-08, 8.490000000000000e-09,9.850000000000000e-09, 1.849000000000000e-08, 6.180000000000000e-09, 4.020000000000000e-09]
	ctor= [3.239000000000000e-08,7.439999999999999e-09,3.090000000000000e-09, 4.719999999999999e-09,1.210000000000000e-09,4.120000000000000e-09, 9.239999999999999e-09,1.027000000000000e-08,6.300000000000000e-09, 8.510000000000000e-09]
	haxi= [1,2,2.900000000000000,3.880000000000000,4.290000000000000,4.850000000000000,5.400000000000000,5.810000000000000,6.160000000000000]
	hrad= [1,2,3,4,4.430000000000000,5,5.380000000000000,5.600000000000000,6, 6.390000000000000]
	htor= [1,2,2.900000000000000,3.880000000000000,4,4.430000000000000, 5.200000000000000,5.400000000000000,5.600000000000000,5.810000000000000]
	zeta1=0.2
	a=0.2
	n=3
	h=[1.0,2.0,3.0,4.0,5.0,6.0]
	C=[
		[42,0.95,4.1,2.75,0.9,1.2],
		[0.2,0.09,0.25,1.0,5.0,0.4]
	]
	Nsurf=10
	D=2.8
	BP=0.4
	PH=1.0e10
	R0=4*60
	mass=[
		[752.094612601836,1.03216046820620e-15,2.23432383705813e-15,3.23726195090639e-17,1.47160061914065e-13,4.28546087505310e-14],
		[5.79397640976254e-16,752.094612601836,-4.27435864480685e-14,-1.04638520070921e-13,2.43945523111095e-15,-6.75674793892966e-15],
		[2.50494069931051e-15,-4.03219124756049e-14,752.094612601836,-1.60062087509864e-13,4.98039109952941e-15,-1.69109226052856e-15],
		[-7.47073394395928e-16,-1.47049039611602e-13,-1.31672450720544e-13,1720.96055133990,-10.5270557086614,13.3850421894362],
		[1.47160061914065e-13,2.62130691004106e-15,4.98039109952941e-15,-10.5270557086614,1352.56559956577,228.801871490305],
		[4.28546087505310e-14,-6.75674793892966e-15,-2.65496910573658e-15,13.3850421894362,228.801871490305,1508.50841532254]
	]
	FgsNom=30.0
	
	#x dict random
	Ru=3000.0
	fc=30.0
	Tst=20.0
	Srg=3e-14
	Sst=2
	Tgs=0.04
	lambda_=1e-6
	Ro=0.98
	QE=0.8
	Mgs=15
	fca=0.01
	Kc=0.0
	Kcf=2000.0
	
	#theta dict
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
	
	#call matlab function
	fncall_return = eng.nexus_lyapunov(
		matlab.double(Ru),
		matlab.double(Us),
		matlab.double(Ud),
		matlab.double(fc),
		matlab.double(Qc),
		matlab.double(Tst),
		matlab.double(Srg),
		matlab.double(Sst),
		matlab.double(Tgs),
		matlab.double(m_SM),
		matlab.double(m_SMhub),
		matlab.double(I_SMhubt),
		matlab.double(I_SMhuba),
		matlab.double(m_RW),
		matlab.double(K_yPM),
		matlab.double(I_xRWA),
		matlab.double(I_yRWA),
		matlab.double(I_RWt),
		matlab.double(I_RWa),
		matlab.double(m_ISO),
		matlab.double(I_ISOa),
		matlab.double(I_ISOt),
		matlab.double(K_yISO),
		matlab.double(K_xISO),
		matlab.double(m_RWAchx),
		matlab.double(I_bus),
		matlab.double(m_bus),
		matlab.double(m_prop),
		matlab.double(I_propt),
		matlab.double(I_propa),
		matlab.double(m_instr),
		matlab.double(I_i1),
		matlab.double(I_i2),
		matlab.double(I_i3),
		matlab.double(A_sptop),
		matlab.double(D_sp),
		matlab.double(t_sp),
		matlab.double(I_ss),
		matlab.double(K_rad1),
		matlab.double(K_rad2),
		matlab.double(K_rISO),
		matlab.double(K_act1),
		matlab.double(K_act2),
		matlab.double(I_iso),
		matlab.double(K_zpet),
		matlab.double(lambda_),
		matlab.double(Ro),
		matlab.double(QE),
		matlab.double(Mgs),
		matlab.double(fca),
		matlab.double(Kc),
		matlab.double(Kcf),
		matlab.double(nray),
		matlab.double(caxi),
		matlab.double(crad),
		matlab.double(ctor),
		matlab.double(haxi),
		matlab.double(hrad),
		matlab.double(htor),
		matlab.double(zeta1),
		matlab.double(a),
		matlab.double(wheel_locs),
		matlab.double(n),
		matlab.double(h),
		matlab.double(C),
		matlab.double(Nsurf),
		matlab.double(D),
		matlab.double(BP),
		matlab.double(PH),
		matlab.double(R0),
		matlab.double(mass),
		matlab.double(FgsNom),
		matlab.double(K_pm1),
		matlab.double(K_pm3),
		matlab.double(K_pm4),
		matlab.double(K_pm5),
		matlab.double(K_pm6),
		matlab.double(K_act_pm2),
		matlab.double(K_act_pm3),
		matlab.double(K_act_pm4),
		matlab.double(K_act_pm5),
		matlab.double(K_act_pm6),
		matlab.double(K_xpet),
		matlab.double(c_RWA),
		matlab.double(c_RWAI),
		matlab.double(c_SM_act),
		matlab.double(c_PM),
		matlab.double(c_PM_act),
		matlab.double(c_petal),
		matlab.double(zeta_sunshield),
		matlab.double(zeta_isolator),
		matlab.double(zeta_solarpanel),
		1 if verbose else 0 #diagnostics
	)
	
	#jitter = fncall_return[1]
	return fncall_return
	
	
if __name__ == '__main__':  
	#This ought to be done at the problem level
	eng = matlab.engine.start_matlab()
	eng.cd(r'nexus_lyapunov', nargout=0)
	
	Q = nexus_lyapunov_system_model([], [], True)
	print(Q)