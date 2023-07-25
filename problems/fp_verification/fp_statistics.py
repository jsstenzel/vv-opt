import sys

sys.path.append('../..')
from problems.functionals import *
	
"""
theta: [0] gain [1] read noise [2] dark current
"""
def fp_hlva(theta, x):
	#define interest params:
	gain = theta[0]
	rn = theta[1]
	dc = theta[2]
	#define parameters:
	tau = x["tau"] #seconds, exposure time

	QoI = math.sqrt(rn**2 + (tau*dc)**2)
	return QoI
	
"""
theta: [0] gain [1] read noise [2] dark current
"""
def fp_qe_hlva(theta, x):
	#define interest params:
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

	QoI = math.sqrt(rn**2 + (tau*dc)**2) / avg_qe
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
	time_gain = testbed_setup + I_gain * (t_gain + t_buffer)
	
	#dark current experiment time
	t_list = []
	for x in range(d_num):
		C = (d_max - math.exp(d_pow)) / (d_num-1)
		t = C * x**d_pow
		t_list.append(t)
	t_list[0] = dc_t0 #clobber; 100ms baseline exposure assumed
	
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