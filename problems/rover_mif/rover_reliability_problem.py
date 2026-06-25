import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
import csv

sys.path.append('../..')
from problems.problem_definition import *
from uq.gsa_plot import *
from approx.kl_divergence import kl_divergence_2gammas

###First, define the fixed def objects
highloadhr_prior = 800
actuator_event_prior = 50
suspension_cycle_prior = 1000
roughkm_prior = 400
nav_km_prior = 1000
nav_event_prior = 100
heater_hr_prior = 200
power_event_prior = 10000
solar_hr_prior = 10000
thermalsensor_event_prior = 100
softwareevent_testnum_prior = 10000
communications_testnum_prior = 2000

#failure rates
theta_defs = [ 
	#MOBILITY/TRACTION
		["wheelwear_per_roughkm", ["gamma_me", [0.006,roughkm_prior]], "continuous"],
		["stall_per_event", ["gamma_me", [0.0015,actuator_event_prior]], "continuous"],
		["harnessfatigue_per_cycle", ["gamma_me", [0.00003,suspension_cycle_prior]], "continuous"],
		["motordrivefault_per_highloadhr", ["gamma_me", [0.0007, highloadhr_prior]], "continuous"],
	#NAVIGATION
		["featureloss_per_lowtexturekm", ["gamma_me", [0.0025,nav_km_prior]], "continuous"],
		["imufault_per_km", ["gamma_me", [0.0018,nav_km_prior]], "continuous"],
		["mapfail_per_attempt", ["gamma_me", [0.003,nav_event_prior]], "continuous"],
		["startrackerfail_per_opportunity", ["gamma_me", [0.001,nav_event_prior]], "continuous"],
		["radiometricdrop_per_opportunity", ["gamma_me", [0.0008,nav_event_prior]], "continuous"],
	#POWER/THERMAL
		["heaterfault_per_coldhr", ["gamma_me", [0.0014,heater_hr_prior]], "continuous"],
		["batteryfault_per_highloadhr", ["gamma_me", [0.0016,highloadhr_prior]], "continuous"],
		["thermalsensorfault_per_opportunity", ["gamma_me", [0.0008,thermalsensor_event_prior]], "continuous"],
		["busfault_per_event", ["gamma_me", [0.001,power_event_prior]], "continuous"],
		["solararrayfault_per_chargehr", ["gamma_me", [0.0005,solar_hr_prior]], "continuous"],
	#COMMUNICATIONS
		["antennafault_per_opportunity", ["gamma_me", [0.0012,communications_testnum_prior]], "continuous"],
		["transceiverfault_per_session", ["gamma_me", [0.0009,communications_testnum_prior]], "continuous"],
		["linkloss_per_window", ["gamma_me", [0.0015,communications_testnum_prior]], "continuous"],
		["grounddropout_per_pass", ["gamma_me", [0.0008,communications_testnum_prior]], "continuous"],
		["commandcorruption_per_event", ["gamma_me", [0.0006,communications_testnum_prior]], "continuous"],
	#SOFTWARE/COMPUTE
		["fswexception_per_cycle", ["gamma_me", [0.0011,softwareevent_testnum_prior]], "continuous"],
		["memcorrupt_per_operation", ["gamma_me", [0.0009,softwareevent_testnum_prior]], "continuous"],
		["watchdogtrip_per_window", ["gamma_me", [0.001,softwareevent_testnum_prior]], "continuous"],
		["invalidplan_per_event", ["gamma_me", [0.0012,softwareevent_testnum_prior]], "continuous"],
		["filterreset_per_update", ["gamma_me", [0.0008,softwareevent_testnum_prior]], "continuous"],
	]

#detected faults
y_defs = [
	#MOBILITY/TRACTION
		"ymeas_wheelwear",
		"ymeas_stall",
		"ymeas_harnessfatigue",
		"ymeas_motordrivefault",
	#NAVIGATION
		"ymeas_featureloss",
		"ymeas_imufault",
		"ymeas_mapfail",
		"ymeas_startrackerfail",
		"ymeas_radiometricdrop",
	#POWER/THERMAL
		"ymeas_heaterfault",
		"ymeas_batteryfault",
		"ymeas_thermalsensorfault",
		"ymeas_busfault",
		"ymeas_solararrayfault",
	#COMMUNICATIONS
		"ymeas_antennafault",
		"ymeas_transceiverfault",
		"ymeas_linkloss",
		"ymeas_grounddropout",
		"ymeas_commandcorruption",
	#SOFTWARE/COMPUTE
		"ymeas_fswexception",
		"ymeas_memcorrupt",
		"ymeas_watchdogtrip",
		"ymeas_invalidplan",
		"ymeas_filterreset",
	]

#rover test drive time for each modal test
d_defs = [
	#MOBILITY/TRACTION
		["test_roughkm_for_wheelwear",				['uniform', [0, 10000]], "continuous"],
		["test_event_for_stall", 					['uniform', [0, 10000]], "discrete"],
		["test_cycle_for_harnessfatigue",			['uniform', [0, 10000]], "discrete"],
		["test_highloadhr_for_motordrivefault", 	['uniform', [0, 10000]], "continuous"],
	#NAVIGATION
		["test_lowtexturekm_for_featureloss",		['uniform', [0, 10000]], "continuous"],
		["test_km_for_imufault",					['uniform', [0, 10000]], "continuous"],
		["test_attempt_for_mapfail", 				['uniform', [0, 10000]], "discrete"],
		["test_opportunity_for_startrackerfail",	['uniform', [0, 10000]], "discrete"],
		["test_opportunity_for_radiometricdrop",	['uniform', [0, 10000]], "discrete"],
	#POWER/THERMAL
		["test_coldhr_for_heaterfault",			['uniform', [0, 10000]], "continuous"],
		["test_highloadhr_for_batteryfault",		['uniform', [0, 10000]], "continuous"],
		["test_opportunity_for_thermalsensorfault",	['uniform', [0, 10000]], "discrete"],
		["test_event_for_busfault",					['uniform', [0, 10000]], "discrete"],
		["test_chargehr_for_solararrayfault", 		['uniform', [0, 10000]], "continuous"],
	#COMMUNICATIONS
		["test_opportunity_for_antennafault",		['uniform', [0, 10000]], "discrete"],
		["test_session_for_transceiverfault",		['uniform', [0, 10000]], "discrete"],
		["test_window_for_linkloss",				['uniform', [0, 10000]], "discrete"],
		["test_pass_for_grounddropout",				['uniform', [0, 10000]], "discrete"],
		["test_event_for_commandcorruption",		['uniform', [0, 10000]], "discrete"],
	#SOFTWARE/COMPUTE
		["test_cycle_for_fswexception",				['uniform', [0, 10000]], "discrete"],
		["test_operation_for_memcorrupt",			['uniform', [0, 10000]], "discrete"],
		["test_window_for_watchdogtrip",			['uniform', [0, 10000]], "discrete"],
		["test_event_for_invalidplan",				['uniform', [0, 10000]], "discrete"],
		["test_update_for_filterreset",				['uniform', [0, 10000]], "discrete"],
	]

#exposure in a 100km segment x interrupt probability x recovery time [hr]
x_defs = [
	#MOBILITY/TRACTION
		["h_wheelwear_per_roughkm", [], "continuous", 				25*0.7*6.0],
		["h_stall_per_event", [], "continuous", 					60*0.054*3.5],
		["h_harnessfatigue_per_cycle", [], "continuous", 			1000*0.5*8.0],
		["h_motordrivefault_per_highloadhr", [], "continuous", 		100*0.75*10.0],
	#NAVIGATION
		["h_featureloss_per_lowtexturekm", [], "continuous", 		30*0.6*4.0],
		["h_imufault_per_km", [], "continuous", 					40*0.7*6.0],
		["h_mapfail_per_attempt", [], "continuous", 				18*0.8*5.0],
		["h_startrackerfail_per_opportunity", [], "continuous", 	12*0.5*2.0],
		["h_radiometricdrop_per_opportunity", [], "continuous", 	10*0.4*8.0],
	#POWER/THERMAL
		["h_heaterfault_per_coldhr", [], "continuous", 				18*0.8*5.0],
		["h_batteryfault_per_highloadhr", [], "continuous", 		14*0.75*6.5],
		["h_thermalsensorfault_per_opportunity", [], "continuous",	22*0.6*4.0],
		["h_busfault_per_event", [], "continuous", 					16*0.85*8.0],
		["h_solararrayfault_per_chargehr", [], "continuous", 		30*0.5*3.0],
	#COMMUNICATIONS
		["h_antennafault_per_opportunity", [], "continuous", 		20*0.7*4.0],
		["h_transceiverfault_per_session", [], "continuous", 		24*0.85*6.0],
		["h_linkloss_per_window", [], "continuous", 				15*0.9*3.5],
		["h_grounddropout_per_pass", [], "continuous", 				10*1.0*8.0],
		["h_commandcorruption_per_event", [], "continuous", 		12*0.95*5.0],
	#SOFTWARE/COMPUTE
		["h_fswexception_per_cycle", [], "continuous", 				28*0.8*5.5],
		["h_memcorrupt_per_operation", [], "continuous", 			18*0.9*7.0],
		["h_watchdogtrip_per_window", [], "continuous", 			16*0.85*6.0],
		["h_invalidplan_per_event", [], "continuous", 				14*0.75*4.5],
		["h_filterreset_per_update", [], "continuous", 				20*0.7*3.5],
	]

###Then, define the possible benchmarking models

#note that y is discrete, result of a Poisson process
#for each theta, grab the d that pairs with it, Poisson -> y result
def rover_drive_test_each(theta, design, x, priormean=None, err=True):
	y = [0]*len(theta)
	for m,(theta_str, tm) in enumerate(theta.items()):
		gen_string = theta_str.split('_')
		dm = design["test_"+gen_string[2]+'_for_'+gen_string[0]]
		
		#sample from Poisson process
		ym = np.random.poisson(tm*dm, 1)[0]
		y[m] = ym
	
	return y

#Q = SUM_m theta_m h_m
def H_total_downtime(theta, x, verbose=False):
	Q = 0
	for m,(theta_str, tm) in enumerate(theta.items()):
		hm = x["h_"+theta_str]
		
		Q += tm*hm
	
	return Q

#add up all time of all tests
def cost_test_time(d, x=None):
	total_time = 0
	
	km_per_hr = 0.5
	software_tests_per_hr = 100
	comms_tests_per_hr = 10
	bench_test_time = 0.25
	mech_test_time = 0.5
	
	#MOBILITY/TRACTION
	total_time += d["test_roughkm_for_wheelwear"] / km_per_hr
	total_time += d["test_event_for_stall"] * mech_test_time
	total_time += d["test_cycle_for_harnessfatigue"] * mech_test_time
	total_time += d["test_highloadhr_for_motordrivefault"] * 1
	#NAVIGATION
	total_time += d["test_lowtexturekm_for_featureloss"] / km_per_hr
	total_time += d["test_km_for_imufault"] / km_per_hr
	total_time += d["test_attempt_for_mapfail"] * mech_test_time
	total_time += d["test_opportunity_for_startrackerfail"] * bench_test_time
	total_time += d["test_opportunity_for_radiometricdrop"] * bench_test_time
	#POWER/THERMAL
	total_time += d["test_coldhr_for_heaterfault"] * 1
	total_time += d["test_highloadhr_for_batteryfault"] * 1
	total_time += d["test_opportunity_for_thermalsensorfault"] * bench_test_time
	total_time += d["test_event_for_busfault"] * bench_test_time
	total_time += d["test_chargehr_for_solararrayfault"] * 1
	#COMMUNICATIONS
	total_time += d["test_opportunity_for_antennafault"] / comms_tests_per_hr
	total_time += d["test_session_for_transceiverfault"] / comms_tests_per_hr
	total_time += d["test_window_for_linkloss"] / comms_tests_per_hr
	total_time += d["test_pass_for_grounddropout"] / comms_tests_per_hr
	total_time += d["test_event_for_commandcorruption"] / comms_tests_per_hr
	#SOFTWARE/COMPUTE
	total_time += d["test_cycle_for_fswexception"] / software_tests_per_hr
	total_time += d["test_operation_for_memcorrupt"] / software_tests_per_hr
	total_time += d["test_window_for_watchdogtrip"] / software_tests_per_hr
	total_time += d["test_event_for_invalidplan"] / software_tests_per_hr
	total_time += d["test_update_for_filterreset"] / software_tests_per_hr
	
	return total_time
	
class ProblemDefinitionPoissonGaussian(ProblemDefinition):
	def __init__(self, _eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs):
		super().__init__(_eta, _H, _G, _theta_defs, _y_defs, _d_defs, _x_defs)
		
		#TODO check to enforce Gamma priors
		
	def sample_posterior(self, y, d):
		vals = [] #a list length num_vals of random numbers of size dim_theta
		for prior in self.priors: ###iterate over dim_theta
			dtype = prior[0]
			params = prior[1]
			#need to do this carefully, we have multiple thetas and multiple samples
	
			if dtype == 'gamma_ab':
				alpha_post = params[0] + y
				beta_post = params[1] + d
				post_thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha_post, scale=1.0/beta_post)
				vals.append(post_thetas_i.tolist())
			elif dtype == 'gamma_mv':
				mean = params[0]
				variance = params[1]
				alpha_post = mean**2 / variance + y
				beta_post = mean / variance + d
				post_thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha_post, scale=1.0/beta_post)
				vals.append(post_thetas_i.tolist())
			elif dtype == 'gamma_me':
				mean = params[0]
				beta = params[1]
				beta_post = beta + d
				#variance = mean / beta
				alpha = mean * beta
				alpha_post = alpha + y
				post_thetas_i = scipy.stats.gamma.rvs(size=num_vals, a=alpha_post, scale=1.0/beta_post)
				vals.append(post_thetas_i.tolist())
			else:
				print("Error, non-Gamma prior")
				sys.exit()

		#turn from list of rvs at each prior, to list of theta rvs
		vals = np.transpose(vals)
	
		#return
		if len(vals) == 1:
			return vals[0] #this is a list (dim_theta)
		else:
			return vals #this is a list (num_vals) of random variables (dim_theta)

	#sum up all of the (weighted) information gains across all modes for a single y
	def infocriterion(self, d, y, rescale_to_downtime=True):
		u=0
		for m,prior in enumerate(self.priors): ###iterate over dim_theta
			dm = d[m]
			ym = y[m]
			dtype = prior[0]
			params = prior[1]
			#need to do this carefully, we have multiple thetas and multiple samples
	
			if dtype == 'gamma_ab':
				am = params[0]
				bm = params[1]
			elif dtype == 'gamma_mv':
				mean = params[0]
				variance = params[1]
				am = mean**2 / variance + y
				bm = mean / variance + d
			elif dtype == 'gamma_me':
				mean = params[0]
				bm = params[1]
				am = mean * bm
			
			#um = am*np.log((bm+dm)/bm) - scipy.special.gammaln(am+ym) + scipy.special.gammaln(am) + ym*scipy.special.digamma(am+ym) - dm*hm*((am+ym)/(bm+dm))
			um = kl_divergence_2gammas(a1=am, b1=bm, a2=am+ym, b2=bm+dm)
			u += um
		return u
		
	def utility(self, d, n_MC):
		thetas = self.prior_rvs(n_MC)
		if n_MC == 1:#ugh edge case
			thetas = [thetas]
		ys = [self.eta(theta, d) for theta in thetas]
		us = [0]*n_MC
		for i,y in enumerate(ys):
			ui = self.infocriterion(d, y) #different utility for each theta,y
			us[i] = ui
		U = np.mean(us)
		return U, us
		
	#Perform an exact analytical sensitivity analysis of H(theta)
	def sensitivity_analysis(self, x=None, doPlot=True):
		if not x:
			x = self.x_default
		S_indices = [None]*self.dim_theta
		labels = self.theta_names
		
		for m, (h_m, theta_m) in enumerate(zip(x, self.priors)):
			#["gamma_me", [0.0011,softwareevent_testnum_prior]]
			dtype = theta_m[0]
			params = theta_m[1]
	
			if dtype == 'gamma_ab':
				alpha = params[0]
				beta = params[1]
				variance = alpha / beta**2
			elif dtype == 'gamma_mv':
				variance = params[1]
			elif dtype == 'gamma_me':
				mean = params[0]
				beta = params[1]
				variance = mean / beta
			else:
				print("Error, non-Gamma prior")
				sys.exit()
			
			s_m = h_m**2 * variance
			S_indices[m] = s_m
		
		s_total = sum(S_indices)
		S_indices = [s_m / s_total for s_m in S_indices]
		
		if doPlot:
			#TODO call the already existing nice plotter
			plot_gsa_full(varnames=labels, S1=S_indices, ST=S_indices, S1_conf=[], ST_conf=[], title="Sensitivity Analysis for System Performance", coplot=False, screening=0, xspin=True)
		
		return S_indices, labels
		

#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
def make_rover_reliability_problem():
	eta = rover_drive_test_each
	Gamma = cost_test_time
	tdefs = theta_defs
	ydefs = y_defs
	ddefs = d_defs
	xdefs = x_defs
	H = H_total_downtime
		
	problem = ProblemDefinitionPoissonGaussian(eta, H, Gamma, tdefs, ydefs, d_defs, x_defs)

	return problem