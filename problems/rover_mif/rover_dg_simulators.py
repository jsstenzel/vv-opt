import sys
import math
import random
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt


sys.path.append('../..')
from models.graphs import *
	
def battery_state_simulator(prev_val):
	#sample from a categorical distribution
	"""
	0 : nominal state
	1 : charging
	2 : high load
	"""
	categories = ["switch_nominal","switch_charging","switch_highload","stay"]
	inertia = (2-prev_val)*0.25 + 1
	weights = [0.6,0.25,0.1,inertia]
	
	sample = random.choices(categories, weights=weights)[0]
	if sample == "switch_nominal":
		return 0
	elif sample == "switch_charging":
		return 1
	elif sample == "switch_highload":
		return 2
	else:
		return prev_val

def _rbf_kernel(X1, X2, sigma_f, l):
        """Computes the RBF covariance between two matrices."""
        dist_matrix = np.sum(X1**2, 1).reshape(-1, 1) + np.sum(X2**2, 1) - 2 * np.dot(X1, X2.T)
        return (sigma_f**2) * np.exp(-0.5 / (l**2) * dist_matrix)

def env_slope_sinemodel(t):
	#Could be modeled with a gaussian process?
	#I'll just model with a sine wave
	#parameters
	max_slope_angle = 10 #degrees
	period = 10 #hours
	
	max_slope = np.tan(np.radians(max_slope_angle))
	amplitude = max_slope / period
	slope = - amplitude * period * np.sin(t / period) #a cosine landscape
	slope_angle = np.tanh(slope)
	
	return slope_angle #radians

def env_sinkage_sampler(t, dt, prev_val, sinkage_mean, sinkage_stddev):
	#Could be modeled as a gaussian process with values forced to be >0?
	#I'll just model with a new Gamma fn sample every hour
	
	"""
	#first, calculate mean sinkage based on physics:
	#parameters
	vehicle_mass = 500 #kg
	wheel_radius = 0.2 #m
	b = 0.1 #m wheel width
	A = b * 0.02 * 4 #rough estimate?
	#properties of lunar surface from #https://www.lpi.usra.edu/publications/books/lunar_sourcebook/pdf/Chapter09.pdf
	n = 1.0 #exponent of soil deformation
	k_c = 1400 #N/m2
	k_phi = 830000 #N/m3 frictional modulus of soil deformation
	k = (k_c/b) + k_phi
	lunar_g = 1.62

	#Bekker's equations (1969)
	W = vehicle_mass * lunar_g
	z_mean = (W / (A*k))**(1/n) #810
	"""
	z_mean = sinkage_mean #0.02 #2 cm, keep it simple
	
	#check if its time (i.e. if we entered a new period in the previous timestep):
	period = 1.0
	if (dt <= period and (t % period) < dt) or dt > period:
		#if timestep is larger than period, then we are always in a new period
		#if so, return a new sample
		sigma = sinkage_stddev #0.02
		M = np.log(z_mean) - sigma**2/2
		V = np.sqrt(np.log(1+sigma**2/z_mean**2))
		return np.random.lognormal(mean=M, sigma=V)
	else:
		#if not, return prev
		return prev_val

def constant_fn(prev_val):
	return prev_val
	
def poisson_process(dt, lambd):
	n_new_events = np.random.poisson(lambd * dt)
	return n_new_events

"""
def env_contamination_PP(dt):
	#parameters
	lambd = .1 #contaminants/hr
	
	#Straightforward Poisson process
	#calculate event count over interval dt
	#return event count at time t+dt
	n_new_events = np.random.poisson(lambd * dt)
	return n_new_events

def env_wheelsnag_PP(dt):
	#parameters
	lambd = 0.01 #snags/hr

	#Straightforward Poisson process
	#calculate event count over interval dt
	#return event count at time t+dt
	n_new_events = np.random.poisson(lambd * dt)
	return n_new_events

def env_rockcontact_PP(dt):
	#parameters
	lambd = 1 #rocks/hr

	#Straightforward Poisson process
	#calculate event count over interval dt
	#return event count at time t+dt
	n_new_events = np.random.poisson(lambd * dt)
	return n_new_events
"""
	
def wheel_stall_NHPP(wheelsnag_hist, rockcontact_hist, wheel_load_hist, dt, lambda_0, snag_dependence, rock_dependence, load_dependence):
	#parameters
	#lambda_0 = 1e-6
	#B1 = 1e-3
	#B2 = 1e-3
	#B3 = 1e-7
	
	N_wheelsnag = sum(wheelsnag_hist)
	N_rockcontact = sum(rockcontact_hist)
	cum_lifetime_wheel_load = sum(wheel_load_hist)
	log_lambd = snag_dependence*N_wheelsnag + rock_dependence*N_rockcontact + load_dependence*cum_lifetime_wheel_load
	n_new_events = np.random.poisson(lambda_0*np.exp(log_lambd-1) * dt)
	return n_new_events

#simulating the torque of the wheel at a time t
def wheel_load_simulator(slope, sinkage, vehicle_mass, wheel_weight_uncertainty, wheel_radius, wheel_width, n, k_c, k_phi):
	#parameters
	"""
	wheel_weight_uncertainty = 10 #N
	vehicle_mass = 500 #kg
	wheel_radius = 0.2 #m
	b = 0.1 #m wheel width
	#properties of lunar surface from #https://www.lpi.usra.edu/publications/books/lunar_sourcebook/pdf/Chapter09.pdf
	n = 1.0 #exponent of soil deformation
	k_c = 1400 #N/m2
	k_phi = 830000 #N/m3 frictional modulus of soil deformation
	"""
	b = wheel_width
	k = (k_c/b) + k_phi
	
	#simulate the physics determined by the instantaneous values for terrain
	#slope - positive and negative slopes determine the load
	lunar_g = 1.62 #m/s2
	weight_on_wheel = (vehicle_mass) * lunar_g * 0.25 * np.random.normal(loc=0, scale=wheel_weight_uncertainty)
	F_slope = weight_on_wheel * math.sin(slope) #slope in radians
	
	#sinkage - nonzero sinkage increases the load
	#Bekker's original equation is for pressure, p = (kc/b + kphi)*z^n https://apps.dtic.mil/sti/tr/pdf/ADA457941.pdf
	#integrate over the depth to get the work
	#multiply by wheel width to get the resistance force:
	F_resistance = (k*b * sinkage**(n+1))/(n+1)
	
	#and then, apply an instantaneous uncertainty to it representing physical uncertainty
	load = max(0,F_slope + F_resistance)
	torque = load * wheel_radius
	
	return torque

def transient_NHPP(battery_state, wheel_load, system_temp, dt, max_nominal_temp, excess_temp_dependence, load_dependence):
	#max_nominal_temp = 273 + 22 #K
	
	#parameters
	#define lambda_0 categorically
	if battery_state == 0:
		lambda_0_battstate = 1e-4
	elif battery_state == 1:
		lambda_0_battstate = 1
	else:
		lambda_0_battstate = 2e-4
	#B1 = 0
	#B2 = 1e-4
	
	excess_temp = max(0, system_temp - max_nominal_temp)
	log_lambd = load_dependence*wheel_load + excess_temp_dependence*excess_temp**2
	n_new_events = np.random.poisson(lambda_0_battstate*np.exp(log_lambd-1) * dt)
	return n_new_events

def brownout_NHPP(battery_state, wheel_load, system_temp, dt, max_nominal_temp, excess_temp_dependence, load_dependence):
	#max_nominal_temp = 273 + 22 #K
	
	#define lambda_0 categorically
	if battery_state == 0:
		lambda_0_battstate = 1e-4
	elif battery_state == 1:
		lambda_0_battstate = 2e-4
	else:
		lambda_0_battstate = 1
	#B1 = 0
	#B2 = 1e-4
	
	excess_temp = max(0, system_temp - max_nominal_temp)
	log_lambd = load_dependence*wheel_load + excess_temp_dependence*excess_temp**2
	n_new_events = np.random.poisson(lambda_0_battstate*np.exp(log_lambd-1) * dt)
	return n_new_events

def system_temperature_simulator(ti_resistance, solar_radiation, battery_state, dt, prev_val, solar_absorptivity, infrared_emissivity, surf_area):
	#parameters
	#solar_absorptivity = 0.15
	#infrared_emissivity = 0.9
	#surf_area = 2.7 * 1.8 #m2
	
	stefan_boltzmann = 5.670374e-8 #W ⋅ m⁻² ⋅ K⁻⁴
	T_prev_K = prev_val #+ 273.14
	
	#calculate the net heat [W]
	Q_solar_radiation = solar_absorptivity * surf_area * solar_radiation
	Q_therm_emission = infrared_emissivity * stefan_boltzmann * surf_area * (T_prev_K**4 - 3**4)
	if battery_state == 0:
		W = 30
	elif battery_state == 1:
		W = 45
	else:
		W = 150
	Q_internal = W
	
	net_heat = Q_solar_radiation + Q_internal - Q_therm_emission
	
	#calculate the dT/dt
	temp_change = net_heat / ti_resistance #this means resistance is K/W
	
	#calculate the temperature at time t
	T = max(0,prev_val + temp_change*dt) #units of Kelvin or deg Celsius
	return T
	

def ti_resistance_simulator(temp_hist, contaminants, prev_val, coeff_thermal_cycling, coeff_contamination):
	#parameters
	#coeff_thermal_cycling = 1e-5 
	#coeff_contamination = 0.0001 #.1% function loss upon each contamination event
	
	#temp_changes = [temp_hist[i+1] - temp_hist[i] for i in range(len(data) - 1)]
	#cum_abs_temp_change = sum([abs(dT) for dT in temp_changes])
	abs_temp_change = abs(temp_hist[-1] - temp_hist[-2]) if len(temp_hist)>1 else 0
	
	#t[0] is nominal, update from prev
	ti_resistance = prev_val 
	#degradation due to thermal cycling
	ti_resistance *= (1 - coeff_thermal_cycling * abs_temp_change)
	#degrade based on contamination events that just happened
	ti_resistance *= (1 - coeff_contamination * contaminants)
	return ti_resistance #units of K/W