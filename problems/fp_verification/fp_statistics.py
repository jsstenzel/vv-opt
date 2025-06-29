import sys
import math

sys.path.append('../..')
from problems.fp_verification.fp_experiment_models import *
from approx.gaussian_process import *
	
"""
theta: [0] gain [1] read noise [2] dark current
"""
def fp_hlva(theta, x, verbose=False):
	#define interest params:
	if type(theta) is dict:
		gain = theta["gain"]
		rn = theta["rn"]
		dc = theta["dc"]
	else:
		gain = theta[0]
		rn = theta[1]
		dc = theta[2]
	#define parameters:
	tau = x["tau"] #seconds, exposure time

	QoI = np.sqrt(rn**2 + (tau*dc))
	return QoI
	
"""
theta: [0] gain [1] read noise [2] dark current
"""
def fp_qe_hlva(theta, x):
	#define interest params:
	if type(theta) is dict:
		gain = theta["gain"]
		rn = theta["rn"]
		dc = theta["dc"]
		qe = theta["qe"]
	else:
		gain = theta[0]
		rn = theta[1]
		dc = theta[2]
		qe = theta[3]
	#define parameters:
	tau = x["tau"] #seconds, exposure time
	
	#we need to turn qe into a single parameter
	#qe is a Functional
	list_wavelength_qe = qe.to_array(10000) #something suitably big
	qes = [wave_qe[1] for wave_qe in list_wavelength_qe]
	avg_qe = np.mean(qes)

	QoI = np.sqrt(rn**2 + (tau*dc)) / avg_qe
	return QoI
	
"""
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
This notional cost model just sums up the times and costs of each isolated decision
"""
def fp_cost_simple(d, x):
	#define design vars:
	t_gain = d["t_gain"]
	I_gain = d["I_gain"]
	n_meas_rn = d["n_meas_rn"]
	d_num = d["d_num"]
	d_max = d["d_max"]
	d_pow = d["d_pow"]
	#define parameters:
	#gain 
	t_gain_setup = x["t_gain_setup"]
	t_gain_buffer = x["t_gain_buffer"]
	#dc
	t_dc_buffer = x["t_dc_buffer"]
	dc_t0 = x["t_0"]
	#rn
	t_rn = x["t_rn"] #exposure time
	t_rn_buffer = x["t_rn_buffer"]
	#cost
	C_engineer = x["C_engineer"]
	testbed_setup = x["testbed_setup"]

	cost = 0
	#gain experiment time
	time_gain = t_gain_setup + (I_gain + 1) * (t_gain + t_gain_buffer) #additional dark exposure on top of I_gain
	
	#dark current experiment time
	t_list = []
	for i in range(d_num):
		t = dark_current_time_fn(i, tmin=dc_t0, dmax=d_max, dpow=d_pow, dnum=d_num)
		t_list.append(t)
	
	time_dc = sum([ti+t_dc_buffer for ti in t_list])

	#read noise experiment time
	time_rn = (t_rn + t_rn_buffer) * n_meas_rn
	
	#time-saving trick: if time_rn = dc_t0, and n_meas_rn <= 3, 
	#you already have the rn data from dc!
	#but if n_meas_rn > 3, you probably want new data all together
	if t_rn == dc_t0 and n_meas_rn <= 3:
		time_rn = 0
	
	#Figure out how these experiment times fit over days
	measurement_time = testbed_setup + time_gain + time_rn + time_dc
	measurement_perday = measurement_time
	measurement_days = 1
	while measurement_perday > 8*3600: #seconds in a workday
		#need to spend another day doing the experiment
		measurement_days += 1
		measurement_time += testbed_setup
		#now see if n days is long enough to do all the experiments
		measurement_perday = measurement_time / measurement_days
	
	cost = C_engineer * measurement_time
	return cost
