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
from problems.rover_mif.rover_dg_simulators import *


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
	"motor winding:[open]" : ["system:(temperature)", "system:[mechanical shock]", "environment:[dust contamination]", "motor winding:(resistance)"],
	"motor winding:[short]" : ["system:(temperature)", "system:[mechanical shock]", "environment:[dust contamination]", "motor winding:(resistance)"],
	"motor winding:(resistance)" : ["system:(temperature)", "system:[mechanical shock]", "environment:[dust contamination]"],
	#DRIVER MOTOR HARNESS FAULT
	"driver-motor harness:[open]" : ["driver-motor harness:(flexure)", "system:(temperature)", "system:(vibration)", "environment:[dust contamination]", "driver-motor harness:(resistance)"],
	"driver-motor harness:[short]" : ["driver-motor harness:(flexure)", "system:(temperature)", "system:(vibration)", "environment:[dust contamination]", "driver-motor harness:(resistance)"],
	"driver-motor harness:(resistance)" : ["driver-motor harness:(flexure)", "system:(temperature)", "system:(vibration)", "environment:[dust contamination]"],
	#THERMAL HARNESS FAULT
	"thermal interface:(resistance)" : ["system:(temperature)", "environment:[dust contamination]"],
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
		["system:(temperature)", "system:[mechanical shock]", "environment:[dust contamination]", "motor winding:(resistance)"])

	#DRIVER MOTOR HARNESS FAULT
	DG_rover.add_nodes_inputs("driver-motor harness:[fault]", #discrete event
		["system:(temperature)", "system:(vibration)", "environment:[dust contamination]"])

	#THERMAL HARNESS FAULT
	DG_rover.add_nodes_inputs("thermal interface:(resistance)", #continuous variable
		["system:(temperature)", "environment:[dust contamination]"])
		
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
		"thermal interface:(resistance)" :["system:(temperature)", "environment:[dust contamination]"], #continuous variable
		#WHEEL FAULT
		"wheel:[stall]" : ["wheel:(load)", "environment:[rock contact]", "environment:[wheel snag]"], #discrete event
		"wheel:(load)" : ["environment:(slope)", "environment:(sinkage)"], #continuous variable
		#POWER BUS FAULT
		"power bus:[transient]" : ["wheel:(load)", "system:(temperature)", "battery:(state)"], #discrete event
		"power bus:[brownout]" : ["wheel:(load)", "system:(temperature)", "battery:(state)"], #discrete event
		#system and environment
		"system:(temperature)" : ["thermal interface:(resistance)", "environment:(solar radiation)", "battery:(state)"] #continuous variable
	}
	DG_rover = SystemDependencyGraph(DG_rover_init, timestep=0.1)
	
	###Set initial conditions
	DG_rover.set_initial_values({
		"environment:(solar radiation)" : 0,
		"thermal interface:(resistance)" : 20000,
		"environment:[dust contamination]" : 0,
		"system:(temperature)" : 230,#10+273, #K
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
	
	###Define constants
	DG_rover.define_constants({
		"CONST_sinkage_mean" : 0.02,
		"CONST_sinkage_stddev" : 0.02, 
		"CONST_contamination_rate" : 0.1, #per hr
		"CONST_rock_rate" : 1, #per hr
		"CONST_snag_rate" : 0.01, #per hr
		"CONST_wheelstall_rate" : 1e-6,
		"CONST_wheelstall_snagdependence" : 1e-3,
		"CONST_wheelstall_rockdependence" : 1e-3,
		"CONST_wheelstall_loaddependence" : 1e-7,
		"CONST_max_nominal_temp" : 273 + 22, #K 
		"CONST_transient_excesstemp_dependence" : 1e-4, 
		"CONST_transient_load_dependence" : 0,
		"CONST_brownout_excesstemp_dependence" : 1e-4, 
		"CONST_brownout_load_dependence" : 0,
		"CONST_solar_absorptivity" : 0.15, 
		"CONST_infrared_emissivity" : 0.9, 
		"CONST_surface_area" : 2.7 * 1.8, #m2
		"CONST_thermresistance_thermalcycling_dependence" : 1e-5,
		"CONST_thermresistance_contamination_dependence" : .0001,
		"CONST_vehicle_mass" : 500, #N 
		"CONST_wheel_weight_sigma" : 10, #N 
		"CONST_wheel_radius" : 0.2, 
		"CONST_wheel_width" : 0.1, 
		"CONST_soil_exponent" : 1.0, 
		"CONST_soil_kc" : 1400, #N/m2
		"CONST_soil_kphi" : 830000, #N/m3 frictional modulus of soil deformation
	})	
	
	
	###Define functions for each node
	#DG_rover.generate_specification_template()
	DG_rover.specify_node("system:(temperature)", system_temperature_simulator, {
        "ti_resistance" : "thermal interface:(resistance)",
        "solar_radiation": "environment:(solar radiation)", 
		"battery_state" : "battery:(state)",
		"solar_absorptivity" : "CONST_solar_absorptivity", 
		"infrared_emissivity" : "CONST_infrared_emissivity", 
		"surf_area" : "CONST_surface_area"
	})
	DG_rover.specify_node("environment:(solar radiation)", constant_fn)
	DG_rover.specify_node("thermal interface:(resistance)", ti_resistance_simulator, 
	{
        "temp_hist" : "system:(temperature)",
        "contaminants" : "environment:[dust contamination]",
		"coeff_thermal_cycling" : "CONST_thermresistance_thermalcycling_dependence",
		"coeff_contamination" : "CONST_thermresistance_contamination_dependence"
	}, get_history=["system:(temperature)"])
	DG_rover.specify_node("environment:[dust contamination]", poisson_process, {
		"lambd" : "CONST_contamination_rate"
	})
	DG_rover.specify_node("wheel:(load)", wheel_load_simulator, {
        "slope" : "environment:(slope)",
        "sinkage" : "environment:(sinkage)",
		"vehicle_mass" : "CONST_vehicle_mass", 
		"wheel_weight_uncertainty" : "CONST_wheel_weight_sigma", 
		"wheel_radius" : "CONST_wheel_radius", 
		"wheel_width" : "CONST_wheel_width", 
		"n" : "CONST_soil_exponent", 
		"k_c" : "CONST_soil_kc", 
		"k_phi" : "CONST_soil_kphi"
	})
	DG_rover.specify_node("wheel:[stall]", wheel_stall_NHPP, 
	{
        "wheel_load_hist" : "wheel:(load)",
        "rockcontact_hist" : "environment:[rock contact]",
        "wheelsnag_hist" : "environment:[wheel snag]",
		"lambda_0" : "CONST_wheelstall_rate",
		"snag_dependence" : "CONST_wheelstall_snagdependence",
		"rock_dependence" : "CONST_wheelstall_rockdependence",
		"load_dependence" : "CONST_wheelstall_loaddependence",
	}, get_history=["wheel:(load)","environment:[rock contact]","environment:[wheel snag]"])
	DG_rover.specify_node("environment:[rock contact]", poisson_process, {
		"lambd" : "CONST_rock_rate"
	})
	DG_rover.specify_node("environment:[wheel snag]", poisson_process, {
		"lambd" : "CONST_snag_rate"
	})
	DG_rover.specify_node("environment:(slope)", env_slope_sinemodel)
	DG_rover.specify_node("environment:(sinkage)", env_sinkage_sampler, {
		"sinkage_mean" : "CONST_sinkage_mean",
		"sinkage_stddev" : "CONST_sinkage_stddev" 
	})
	DG_rover.specify_node("power bus:[transient]", transient_NHPP, {
        "wheel_load" : "wheel:(load)",
        "system_temp" : "system:(temperature)",
        "battery_state" : "battery:(state)",
		"max_nominal_temp" : "CONST_max_nominal_temp", 
		"excess_temp_dependence" : "CONST_transient_excesstemp_dependence", 
		"load_dependence" : "CONST_transient_load_dependence"
	})
	DG_rover.specify_node("battery:(state)", battery_state_simulator)
	DG_rover.specify_node("power bus:[brownout]", brownout_NHPP, {
        "wheel_load" : "wheel:(load)",
        "system_temp" : "system:(temperature)",
        "battery_state" : "battery:(state)",
		"max_nominal_temp" : "CONST_max_nominal_temp", 
		"excess_temp_dependence" : "CONST_brownout_excesstemp_dependence", 
		"load_dependence" : "CONST_brownout_load_dependence"
	})

	#testing
	#DG_rover.display_interactive()
	
	DG_rover.get_specifications()
	#DG_rover.printout()
	
	"""
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
	"""
	times,values_timeseries,last_values = DG_rover.run_simulation(10, startOver=True, doPrint=True, logfile="rover_data")
		
	node_value_t = np.array(values_timeseries).T.tolist()
	fig, axes = plt.subplots(nrows=len(node_value_t), ncols=1, figsize=(6, 2 * len(node_value_t)), sharex=True)
	for i, ax in enumerate(axes):
		ax.plot(times, node_value_t[i])
		ax.set_ylabel(
			last_values[i][0],
			rotation=0,       # Forces text to be horizontal
			labelpad=20,      # Adds spacing so text doesn't overlap tick marks
			va='center',      # Centers text vertically relative to the axis
			ha='right'        # Aligns text cleanly to the right edge
		)
	plt.xlabel('t [hours]')
	plt.tight_layout()  # Prevents overlapping labels and titles
	plt.show()
	

	
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