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



def full_problem():
	#process:
	#- identify the events/parameters in each failure mode
	#- then in each cause/driver
	#- then in each observable symptom & local effect

	#note:
	#- all variables have a time-dependence in general; think of them all as time series
	#- in addition to this interdependence, discrete events can also have an intrinsic constant or time-varying failure rate

	DG_rover_init = {
	#POWER INVERTER FAULT
	"power inverter:[overtemp fault]" : ["system:(current)", "system:(duty cycle)", "system:(temperature)"],
	#POWER ELECTRONICS SWITCHING DEVICE FAULT
	"power electronics switching device:[open]" : ["system:(temperature)", "system:(current)", "system:[electrical transient]", "environment:(radiation)"],
	"power electronics switching device:[short]" : ["system:(temperature)", "system:(current)", "system:[electrical transient]", "environment:(radiation)"],
	"power electronics switching device:(degradation)" : ["system:(temperature)", "system:(current)", "system:[electrical transient]", "environment:(radiation)"],
	#CURRENT SENSING PATH FAULT
	"current sensing path:[open]" : ["system:(temperature)"],
	"current sensing path:(bias)" : ["ADC:(drift)", "system:(temperature)"],
	"current sensing path:(drift)" : ["ADC:(drift)", "system:(temperature)"],
	"current sensing path:(saturation)" : ["ADC:(drift)", "system:(temperature)"],
	#HALL SENSOR FAULT
	"hall sensor:[invalid signal fault]" : ["system:(temperature)", "environment:(radiation)"],
	#FPGA FAULT
	"motor controller/FPGA:[timing fault]" : ["wheel:(load)"],
	"motor controller/FPGA:(saturation)" : ["wheel:(load)"],
	#MOTOR WINDING FAULT
	"motor winding:[open]" : ["system:(temperature)", "system:[mechanical shock]", "environment:[contamination]", "motor winding:(resistance)"],
	"motor winding:[short]" : ["system:(temperature)", "system:[mechanical shock]", "environment:[contamination]", "motor winding:(resistance)"],
	"motor winding:(resistance)" : ["system:(temperature)", "system:[mechanical shock]", "environment:[contamination]"],
	#DRIVER MOTOR HARNESS FAULT
	"driver-motor harness:[open]" : ["driver-motor harness:(flexure)", "system:(temperature)", "system:(vibration)", "environment:[contamination]", "driver-motor harness:(resistance)"],
	"driver-motor harness:[short]" : ["driver-motor harness:(flexure)", "system:(temperature)", "system:(vibration)", "environment:[contamination]", "driver-motor harness:(resistance)"],
	"driver-motor harness:(resistance)" : ["driver-motor harness:(flexure)", "system:(temperature)", "system:(vibration)", "environment:[contamination]"],
	#THERMAL HARNESS FAULT
	"thermal interface:(resistance)" : ["system:(temperature)", "environment:[contamination]"],
	"system:(temperature)" : ["thermal interface:(resistance)", "environment:(temperature)"],
	#WHEEL FAULT
	"wheel:[stall]" : ["wheel:(load)", "environment:[rock contact]", "environment:[wheel snag]"],
	"wheel:(load)" : ["environment:(slope)", "environment:(sinkage)"],
	#POWER BUS FAULT
	"power bus:[transient]" : ["wheel:(load)", "system:(temperature)", "battery:(state)"],
	"power bus:[brownout]" : ["wheel:(load)", "system:(temperature)", "battery:(state)"],
	#ALSO THINK ABOUT THE ENVIRONMENT:
	#NOW THINK ABOUT WHAT GETS SENSED:
	}
	
	DG_rover = SystemDependencyGraph(DG_rover_init)

	#testing
	#DG_rover.printout()
	#DG_rover.display_simple()
	#DG_rover.display("spring")
	DG_rover.display_interactive()

"""
def simplified_problem():
	DG_rover = SystemDependencyGraph()

	DG_rover.add_nodes_inputs("power inverter:[fault]", #discrete event
		["system:(current)", "system:(duty cycle)", "system:(temperature)"])

	#POWER ELECTRONICS SWITCHING DEVICE FAULT
	DG_rover.add_nodes_inputs("power electronics switching device:[fault]", #discrete event
		["system:(temperature)", "system:(current)", "system:[electrical transient]", "environment:(radiation)"])

	#CURRENT SENSING PATH FAULT
	DG_rover.add_nodes_inputs("current sensing path:[fault]", #continuous variable
		["ADC:(drift)", "system:(temperature)"])

	#HALL SENSOR FAULT
	DG_rover.add_nodes_inputs("hall sensor:[fault]", #discrete event
		["system:(temperature)", "environment:(radiation)"])

	#FPGA FAULT
	DG_rover.add_nodes_inputs("motor controller/FPGA:[fault]", #discrete event
		["wheel:(load)"])

	#MOTOR WINDING FAULT
	DG_rover.add_nodes_inputs("motor winding:[fault]", #discrete event
		["system:(temperature)", "system:[mechanical shock]", "environment:[contamination]", "motor winding:(resistance)"])

	#DRIVER MOTOR HARNESS FAULT
	DG_rover.add_nodes_inputs("driver-motor harness:[fault]", #discrete event
		["system:(temperature)", "system:(vibration)", "environment:[contamination]"])

	#THERMAL HARNESS FAULT
	DG_rover.add_nodes_inputs("thermal interface:(resistance)", #continuous variable
		["system:(temperature)", "environment:[contamination]"])
		
	DG_rover.add_nodes_inputs("system:(temperature)", #continuous variable
		["thermal interface:(resistance)", "environment:(temperature)"])

	#WHEEL FAULT
	DG_rover.add_nodes_inputs("wheel:[fault]", #discrete event
		["environment:[rock contact]", "environment:[wheel snag]", "environment:(slope)", "environment:(sinkage)"])

	#POWER BUS FAULT
	DG_rover.add_nodes_inputs("power bus:[fault]", #discrete event
		["wheel:(load)", "system:(temperature)", "battery:(state)"])

	#testing
	DG_rover.printout()
	DG_rover.display_simple()
	DG_rover.display()
	"""

def subproblem():
	###Define graph structure
	DG_rover_init = {
		#THERMAL HARNESS FAULT
		"thermal interface:(resistance)" :["system:(temperature)", "environment:[contamination]"], #continuous variable
		#WHEEL FAULT
		"wheel:[stall]" : ["wheel:(load)", "environment:[rock contact]", "environment:[wheel snag]"], #discrete event
		"wheel:(load)" : ["environment:(slope)", "environment:(sinkage)"], #continuous variable
		#POWER BUS FAULT
		"power bus:[transient]" : ["wheel:(load)", "system:(temperature)", "battery:(state)"], #discrete event
		"power bus:[brownout]" : ["wheel:(load)", "system:(temperature)", "battery:(state)"], #discrete event
		#system and environment
		"system:(temperature)" : ["thermal interface:(resistance)", "env:(solar radiation)", "battery:(state)"] #continuous variable
	}
	DG_rover = SystemDependencyGraph(DG_rover_init, timestep=0.1)
	
	###Define functions for each node
	#DG_rover.generate_specification_template()
	DG_rover.specify_node("system:(temperature)", system_temperature_simulator, {
        "ti_resistance" : "thermal interface:(resistance)",
        "solar_radiation": "env:(solar radiation)", 
		"battery_state" : "battery:(state)"
	})
	DG_rover.specify_node("env:(solar radiation)", constant_fn)
	DG_rover.specify_node("thermal interface:(resistance)", ti_resistance_simulator, 
	{
        "temp_hist" : "system:(temperature)",
        "contaminants" : "environment:[contamination]",
	}, get_history=["system:(temperature)"])
	DG_rover.specify_node("environment:[contamination]", env_contamination_PP)
	DG_rover.specify_node("wheel:(load)", wheel_load_simulator, {
        "slope" : "environment:(slope)",
        "sinkage" : "environment:(sinkage)",
	})
	DG_rover.specify_node("wheel:[stall]", wheel_stall_NHPP, 
	{
        "wheel_load_hist" : "wheel:(load)",
        "rockcontact_hist" : "environment:[rock contact]",
        "wheelsnag_hist" : "environment:[wheel snag]",
	}, get_history=["wheel:(load)","environment:[rock contact]","environment:[wheel snag]"])
	DG_rover.specify_node("environment:[rock contact]", env_rockcontact_PP)
	DG_rover.specify_node("environment:[wheel snag]", env_wheelsnag_PP)
	DG_rover.specify_node("environment:(slope)", env_slope_GP)
	DG_rover.specify_node("environment:(sinkage)", env_sinkage_GP)
	DG_rover.specify_node("power bus:[transient]", transient_NHPP, {
        "wheel_load" : "wheel:(load)",
        "system_temp" : "system:(temperature)",
        "battery_state" : "battery:(state)",
	})
	DG_rover.specify_node("battery:(state)", battery_state_simulator)
	DG_rover.specify_node("power bus:[brownout]", brownout_NHPP, {
        "wheel_load" : "wheel:(load)",
        "system_temp" : "system:(temperature)",
        "battery_state" : "battery:(state)",
	})
	
	###Set initial conditions
	DG_rover.set_initial_values({
		"env:(solar radiation)" : 0,
		"thermal interface:(resistance)" : 20000,
		"environment:[contamination]" : 0,
		"system:(temperature)" : 22+273, #K
		"wheel:(load)" : 0,
		"wheel:[stall]" : 0,
		"environment:[rock contact]" : 0,
		"environment:[wheel snag]" : 0,
		"environment:(slope)" : 0,
		"environment:(sinkage)" : 0.02, #2 cm
		"power bus:[transient]" : 0,
		"battery:(state)" : 0,
		"power bus:[brownout]" : 0
	})

	#testing
	#DG_rover.display_interactive()
	
	DG_rover.get_specifications()
	#DG_rover.printout()
	
	print("t=0")
	DG_rover.init_simulation()
	DG_rover.get_node_vals(doPrint=True)
	
	#run for X time
	run_time = 72 #hours
	time = 0
	values_timeseries = []
	times=[]
	while time < run_time:
		time, values = DG_rover.simulation_step()
		print("t="+str(time))
		times.append(time)
		values_timeseries.append([value[1] for value in values])
	
	for item in values:
		print(item)
		
	node_value_t = np.array(values_timeseries).T.tolist()
	fig, axes = plt.subplots(nrows=len(node_value_t), ncols=1, figsize=(6, 2 * len(node_value_t)), sharex=True)
	for i, ax in enumerate(axes):
		ax.plot(times, node_value_t[i])
		ax.set_ylabel(
			values[i][0],
			rotation=0,       # Forces text to be horizontal
			labelpad=20,      # Adds spacing so text doesn't overlap tick marks
			va='center',      # Centers text vertically relative to the axis
			ha='right'        # Aligns text cleanly to the right edge
		)
	plt.xlabel('t [hours]')
	plt.tight_layout()  # Prevents overlapping labels and titles
	plt.show()
	
def battery_state_simulator(prev_val):
	#sample from a categorical distribution
	"""
	0 : nominal state
	1 : charging
	2 : high load
	"""
	categories = ["switch_nominal","switch_charging","switch_highload","stay"]
	inertia = (2-prev_val)*0.25
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

def env_slope_GP(t):
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

def env_sinkage_GP(t, dt, prev_val):
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
	z_mean = 0.02 #2 cm, keep it simple
	
	#check if its time (i.e. if we entered a new period in the previous timestep):
	period = 1.0
	if (dt <= period and (t % period) < dt) or dt > period:
		#if timestep is larger than period, then we are always in a new period
		#if so, return a new sample
		sigma = 0.02
		M = np.log(z_mean) - sigma**2/2
		V = np.sqrt(np.log(1+sigma**2/z_mean**2))
		return np.random.lognormal(mean=M, sigma=V)
	else:
		#if not, return prev
		return prev_val

def constant_fn(prev_val):
	return prev_val

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
	lambd = 0.0001 #snags/hr

	#Straightforward Poisson process
	#calculate event count over interval dt
	#return event count at time t+dt
	n_new_events = np.random.poisson(lambd * dt)
	return n_new_events

def env_rockcontact_PP(dt):
	#parameters
	lambd = 0.1 #rocks/hr

	#Straightforward Poisson process
	#calculate event count over interval dt
	#return event count at time t+dt
	n_new_events = np.random.poisson(lambd * dt)
	return n_new_events
	
def wheel_stall_NHPP(wheelsnag_hist, rockcontact_hist, wheel_load_hist, dt):
	#parameters
	lambda_0 = 1e-6
	B1 = 1e-6
	B2 = 1e-6
	B3 = 1e-10
	
	N_wheelsnag = sum(wheelsnag_hist)
	N_rockcontact = sum(rockcontact_hist)
	cum_lifetime_wheel_load = sum(wheel_load_hist)
	log_lambd = B1*N_wheelsnag + B2*N_rockcontact + B3*cum_lifetime_wheel_load
	n_new_events = np.random.poisson(lambda_0*np.exp(log_lambd-1) * dt)
	return n_new_events

#simulating the torque of the wheel at a time t
def wheel_load_simulator(slope, sinkage):
	#parameters
	wheel_weight_uncertainty = 10 #N
	vehicle_mass = 500 #kg
	wheel_radius = 0.2 #m
	b = 0.1 #m wheel width
	#properties of lunar surface from #https://www.lpi.usra.edu/publications/books/lunar_sourcebook/pdf/Chapter09.pdf
	n = 1.0 #exponent of soil deformation
	k_c = 1400 #N/m2
	k_phi = 830000 #N/m3 frictional modulus of soil deformation
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

def transient_NHPP(battery_state, wheel_load, system_temp, dt):
	max_nominal_temp = 273 + 22 #K
	
	#parameters
	#define lambda_0 categorically
	if battery_state == 0:
		lambda_0_battstate = 1e-4
	elif battery_state == 1:
		lambda_0_battstate = 1
	else:
		lambda_0_battstate = 2e-4
	B1 = 0
	B2 = 1e-4
	
	excess_temp = max(0, system_temp - max_nominal_temp)
	log_lambd = B1*wheel_load + B2*excess_temp**2
	n_new_events = np.random.poisson(lambda_0_battstate*np.exp(log_lambd-1) * dt)
	return n_new_events

def brownout_NHPP(battery_state, wheel_load, system_temp, dt):
	max_nominal_temp = 273 + 22 #K
	
	#define lambda_0 categorically
	if battery_state == 0:
		lambda_0_battstate = 1e-4
	elif battery_state == 1:
		lambda_0_battstate = 2e-4
	else:
		lambda_0_battstate = 1
	B1 = 0
	B2 = 1e-4
	
	excess_temp = max(0, system_temp - max_nominal_temp)
	log_lambd =  B1*wheel_load + B2*excess_temp**2
	n_new_events = np.random.poisson(lambda_0_battstate*np.exp(log_lambd-1) * dt)
	return n_new_events

def system_temperature_simulator(ti_resistance, solar_radiation, battery_state, dt, prev_val):
	#parameters
	solar_absorptivity = 0.15
	infrared_emissivity = 0.9
	surf_area = 2.7 * 1.8 #m2
	
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
	

def ti_resistance_simulator(temp_hist, contaminants, prev_val):
	#parameters
	coeff_thermal_cycling = 1e-5 
	coeff_contamination = 0.0001 #.1% function loss upon each contamination event
	
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
	
	

	
if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Function to run for this vvopt analysis')
	#parser.add_argument('-n', type=int, default=0, help='Number of iterations to give to the function')
	#parser.add_argument('-v', type=float, default=0, help='Function value input')
	args = parser.parse_args()
	
	if args.run == "full_problem":
		full_problem()
		
	if args.run == "simplified_problem":
		simplified_problem()
		
	if args.run == "subproblem":
		subproblem()
		
	if args.run == 'debug':
		debug()