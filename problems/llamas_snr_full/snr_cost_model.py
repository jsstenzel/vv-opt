import sys
import math

sys.path.append('../..')
from approx.gaussian_process import *
from problems.llamas_snr_full.snr_exp_models import dark_current_time_fn

	
"""
theta: 
x parameters: 
This notional cost model just sums up the times and costs of each isolated decision
"""
def snr_cost(d, x):
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
	#qe
	#Handle measurement number: you can do 0, or 3 or more
	t_qe = d["t_qe"] if d["t_qe"]>0.1 else 0.1 #times smaller than this just aren't meaningful
	n_qe = d["n_qe"] if (d["n_qe"]==0 or d["n_qe"]>3) else 3
	qe_setup = x["t_qe_setup"]
	qe_buffer = x["t_qe_buffer"]
	#VPH
	#Handle measurement number: you can do 0, or 3 or more
	d_vph_n_pts = d["d_vph_n_pts"] if (d["d_vph_n_pts"]==0 or d["d_vph_n_pts"]>3) else 3
	vph_setup = x["t_vph_setup"]
	vph_exposure_per_pt = x["t_vph_per_pt"]
	#Dichroic
	#Handle measurement number: you can do 0, or 4 or more
	d_dichroic_n_pts = d["d_dichroic_n_pts"] if (d["d_dichroic_n_pts"]==0 or d["d_dichroic_n_pts"]>4) else 4
	dichroic_setup = x["t_dichroic_setup"]
	dichroic_exposure_per_pt = x["t_dichroic_per_pt"]
	#Collimator
	d_coll_n_pts = d["d_coll_n_pts"]
	coll_setup = x["t_coll_setup"]
	coll_exposure_per_pt = x["t_coll_per_pt"]
	#Lens
	d_redcam_n_pts = d["d_redcam_n_pts"]
	d_greencam_n_pts = d["d_greencam_n_pts"]
	d_bluecam_n_pts = d["d_bluecam_n_pts"]
	camera_test_setup = x["t_camera_test_setup"]
	camera_exposure_per_pt = x["t_camera_per_pt"]
	#Fiber
	d_frd_n_meas = d["d_frd_n_meas"]
	frd_setup = x["t_frd_setup"]
	frd_test_time = x["t_frd_test"]
	#cost
	C_engineer = x["C_engineer"]
	fp_testbed_setup = x["fp_testbed_setup"]
	
	cost = 0

	###################################################
	# CCD Experiments
	###################################################
	#gain experiment time
	time_gain = t_gain_setup + (I_gain + 1) * (t_gain + t_gain_buffer) if I_gain>0 else 0.0 #additional dark exposure on top of I_gain
	
	#dark current experiment time
	t_list = []
	if d_num > 1:
		for i in range(d_num):
			t = dark_current_time_fn(i, tmin=dc_t0, dmax=d_max, dpow=d_pow, dnum=d_num)
			t_list.append(t)
	
		time_dc = sum([ti+t_dc_buffer for ti in t_list])
	else:
		time_dc = 0.0

	#read noise experiment time
	time_rn = (t_rn + t_rn_buffer) * n_meas_rn
	
	#time-saving trick: if time_rn = dc_t0, and n_meas_rn <= 3, 
	#you already have the rn data from dc!
	#but if n_meas_rn > 3, you probably want new data all together
	if t_rn == dc_t0 and n_meas_rn <= 3:
		time_rn = 0
		
	###########
	# QE experiment
	time_qe = qe_setup + n_qe * (t_qe + qe_buffer) if n_qe > 0 else 0.0
	
	#Figure out how these experiment times fit over days
	fp_time = time_gain + time_rn + time_dc + time_qe
	if fp_time > 0.0:
		fp_time += fp_testbed_setup #only add setup if you actually do fp experiments!
	measurement_perday = fp_time
	measurement_days = 1
	while measurement_perday > 8*3600: #seconds in a workday
		#need to spend another day doing the experiment
		measurement_days += 1
		fp_time += fp_testbed_setup
		#now see if n days is long enough to do all the experiments
		measurement_perday = fp_time / measurement_days
	
	fp_time *= 3 #for each of the 3 cameras, which require a new setup each time
	
	###################################################
	# Component Measurements
	###################################################
	#All of these are straightforward throughput measurements, except the FRD
	
	#VPH - about a month was scheduled for it in the WBS, Jan 13 2021 to Feb 19 2021
	#Thats 19 work days, or 152 work hours
	#But actually doing it is probably a fraction of that time. Maybe 1/10th even.
	#For 24 VPH and 3 measurements each, over lets say 15 hours, thats... hard to extrapolate from
	
	#VPH - n_pts corresponds to individual measurements
	time_vph = 3 * (vph_setup + d_vph_n_pts * (vph_exposure_per_pt)) if d_vph_n_pts>0 else 0.0
	#Camera - n_pts corresponds to individual measurements
	time_lenses = (camera_test_setup + d_redcam_n_pts * (camera_exposure_per_pt)) if d_redcam_n_pts>0 else 0.0
	time_lenses += (camera_test_setup + d_greencam_n_pts * (camera_exposure_per_pt)) if d_greencam_n_pts>0 else 0.0
	time_lenses += (camera_test_setup + d_bluecam_n_pts * (camera_exposure_per_pt)) if d_bluecam_n_pts>0 else 0.0
	
	#Dichroic - n_pts corresponds to spectral resolution, which trades off with integration time
	#i.e. this all happens in one exposure. that exposure probably has a minimum time of like 1 second
	time_dichroic = 2 * (dichroic_setup + max(d_dichroic_n_pts * dichroic_exposure_per_pt,1)) if d_dichroic_n_pts>0 else 0.0
	#Collimator - n_pts corresponds to spectral resolution, which trades off with integration time
	time_coll = coll_setup + max(d_coll_n_pts * coll_exposure_per_pt,1) if d_coll_n_pts>0 else 0.0

	#Fiber frd - n_pts corresponds to individual measurements
	time_frd = frd_setup + d_frd_n_meas * frd_test_time if d_frd_n_meas>0 else 0.0
	
	cost = C_engineer * (fp_time + time_vph + time_dichroic + time_coll + time_lenses + time_frd) 
	return cost
