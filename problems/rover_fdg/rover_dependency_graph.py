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
from problems.rover_fdg.rover_dg_simulators import *



def extended_problem():
	###Define graph structure
	DG_rover_init = {
		#MISSION
		"mission:(downtime)" : ["control:(fault recovery)"],
		"mission:(distance traversed)": ["mission:(productive traversal)", "mission:(unproductive traversal)"], #simple addition
		"mission:(productive traversal)" : ["control:(operational mode)", "control:(drive speed)", "mission:(unproductive traversal remaining)"], #drive at speed if in drive mode
		#FAULTS
		"mission:(unproductive traversal remaining)" : ["control:[avoid obstacle]", "software:[nav error]"],
		"power bus:[transient]" : ["wheel:(torque)", "WEB:(temperature)"], #discrete event
		"power bus:[brownout]" : ["WEB:(temperature)", "battery:(power draw)"], #discrete event
		"actuator:[stall]" : ["wheel:(torque)", "motor drive electronics:(stress)"], #discrete event, when torque hits a limit
		"actuator:[loss]" : ["wheel:(torque)", "environment:[dust contamination]"],
		"wheel:[snag]" : ["environment:[rock contact]"],
		"motor drive electronics:[fault]" : ["motor drive electronics:(stress)"],
		#FDIR??
		"fdir:[transient detected]" : ["power bus:[transient]"],
		"fdir:[brownout detected]" : ["power bus:[brownout]"],
		"fdir:[stall detected]" : ["actuator:[stall]"],
		"fdir:[actuator loss detected]" : ["actuator:[loss]"],
		"fdir:[mde fault detected]" : ["motor drive electronics:[fault]"],
		"control:(fault recovery)" : ["fdir:[transient detected]", "fdir:[brownout detected]", "fdir:[stall detected]", "fdir:[actuator loss detected]", "fdir:[mde fault detected]"],
		#SYSTEM & CONOPS
		"system:(pose)" : ["mission:(productive traversal)", "mission:(unproductive traversal)"],
		"control:(operational mode)" : ["environment:(UTC time)", "control:(fault recovery)"], #if we're accumulating downtime, then we have to be in safe mode
		"control:(drive speed)" : ["actuator:[loss]", "fdir:[actuator loss detected]"], #with each actuator down, we have to go slower
		"control:[avoid obstacle]" : ["environment:[obstacle]"],
		#ENVIRONMENTAL INTERACTION
		"environment:[dust contamination]" : ["environment:(temperature)", "environment:(slip ratio)"],
		"environment:(temperature)" : ["environment:(UTC time)"],
		"environment:(slip ratio)" : ["system:(pose)", "environment:(slope)", "environment:(sinkage)"],
		"environment:(slope)" : ["system:(pose)"],
		"environment:(sinkage)" : ["system:(pose)", "environment:(slope)"],
		#THERMAL
		"radiator:(efficiency)" :["environment:[dust contamination]"], #continuous variable
		"thermal interface:(resistance)" :["WEB:(temperature)", "environment:[dust contamination]"], #continuous variable
		"WEB:(temperature)" : ["thermal interface:(resistance)", "environment:(temperature)", "battery:(power draw)", "radiator:(efficiency)"], #continuous variable
		#MOBILITY
		"wheel:(torque)" : ["environment:(slope)", "environment:(slip ratio)", "environment:(sinkage)", "environment:[rock contact]", "wheel:[snag]", "wheel:(degradation)", "control:(drive speed)", "actuator:[loss]"], #continuous variable
		"wheel:(degradation)" : ["mission:(distance traversed)", "environment:[rock contact]", "wheel:[snag]"],
		"motor drive electronics:(stress)" : ["wheel:(torque)"],
		#POWER
		"battery:(power draw)" : ["control:(operational mode)"],
	}
	"""
		#POWER INVERTER FAULT
		"power inverter:[overtemp fault]" : ["system:(current)", "system:(duty cycle)", "WEB:(temperature)"],
		#POWER ELECTRONICS SWITCHING DEVICE FAULT
		"power electronics switching device:[open]" : ["WEB:(temperature)", "system:(current)", "system:[electrical transient]", "environment:(radiation)"],
		"power electronics switching device:[short]" : ["WEB:(temperature)", "system:(current)", "system:[electrical transient]", "environment:(radiation)"],
		"power electronics switching device:(degradation)" : ["WEB:(temperature)", "system:(current)", "system:[electrical transient]", "environment:(radiation)"],
		#CURRENT SENSING PATH FAULT
		"current sensing path:[open]" : ["WEB:(temperature)"],
		"current sensing path:(bias)" : ["ADC:(drift)", "WEB:(temperature)"],
		"current sensing path:(drift)" : ["ADC:(drift)", "WEB:(temperature)"],
		"current sensing path:(saturation)" : ["ADC:(drift)", "WEB:(temperature)"],
		#HALL SENSOR FAULT
		"hall sensor:[invalid signal fault]" : ["WEB:(temperature)", "environment:(radiation)"],
		#FPGA FAULT
		"motor controller/FPGA:[timing fault]" : ["wheel:(load)"],
		"motor controller/FPGA:(saturation)" : ["wheel:(load)"],
		#MOTOR WINDING FAULT
		"motor winding:[open]" : ["WEB:(temperature)", "system:[mechanical shock]", "environment:[dust contamination]", "motor winding:(resistance)"],
		"motor winding:[short]" : ["WEB:(temperature)", "system:[mechanical shock]", "environment:[dust contamination]", "motor winding:(resistance)"],
		"motor winding:(resistance)" : ["WEB:(temperature)", "system:[mechanical shock]", "environment:[dust contamination]"],
		#DRIVER MOTOR HARNESS FAULT
		"driver-motor harness:[open]" : ["driver-motor harness:(flexure)", "WEB:(temperature)", "system:(vibration)", "environment:[dust contamination]", "driver-motor harness:(resistance)"],
		"driver-motor harness:[short]" : ["driver-motor harness:(flexure)", "WEB:(temperature)", "system:(vibration)", "environment:[dust contamination]", "driver-motor harness:(resistance)"],
		"driver-motor harness:(resistance)" : ["driver-motor harness:(flexure)", "WEB:(temperature)", "system:(vibration)", "environment:[dust contamination]"],
		"""
	DG_rover = SystemDependencyGraph(DG_rover_init, timestep=0.1)
	DG_rover.printout()
	
	###Set categories and display
	DG_rover_categorization = {
		"control:(fault recovery)" : "CONTROL",
		"mission:(downtime)" : "MISSION",
		"mission:(productive traversal)" : "MISSION",
		"mission:(distance traversed)" : "MISSION",
		"mission:(unproductive traversal)" : "MISSION",
		"control:(operational mode)" : "CONTROL",
		"control:(drive speed)" : "CONTROL",
		"mission:(unproductive traversal)" : "MISSION",
		"control:[avoid obstacle]" : "CONTROL",
		"wheel:(torque)" : "SYSTEM",
		"power bus:[transient]" : "FAULT",
		"WEB:(temperature)" : "SYSTEM",
		"power bus:[brownout]" : "FAULT",
		"battery:(power draw)" : "SYSTEM",
		"actuator:[stall]" : "FAULT",
		"motor drive electronics:(stress)" : "SYSTEM",
		"actuator:[loss]" : "FAULT",
		"environment:[dust contamination]" : "ENVIRONMENT",
		"environment:[rock contact]" : "ENVIRONMENT",
		"wheel:[snag]" : "FAULT",
		"motor drive electronics:[fault]" : "FAULT",
		"fdir:[transient detected]" : "OBSERVATION",
		"fdir:[brownout detected]" : "OBSERVATION",
		"fdir:[stall detected]" : "OBSERVATION",
		"fdir:[actuator loss detected]" : "OBSERVATION",
		"fdir:[mde fault detected]" : "OBSERVATION",
		"system:(pose)" : "MISSION",
		"environment:(UTC time)" : "ENVIRONMENT",
		"environment:[obstacle]" : "FAULT",
		"environment:(temperature)" : "ENVIRONMENT",
		"environment:(slip ratio)" : "ENVIRONMENT",
		"environment:(slope)" : "ENVIRONMENT",
		"environment:(sinkage)" : "ENVIRONMENT",
		"radiator:(efficiency)" : "SYSTEM",
		"thermal interface:(resistance)" : "SYSTEM",
		"wheel:(degradation)" : "SYSTEM",
	}
	DG_rover.categorize_nodes(DG_rover_categorization)
	
	category_colors = {
		"FAULT":"red",
		"MISSION":"black",
		"ENVIRONMENT":"green",
		"SYSTEM":"orange",
		"OBSERVATION":"blue",
		"CONTROL":"purple",
	}
	DG_rover.set_category_colors(category_colors)
	
	DG_rover.display_interactive()
	
	
def restricted_problem():
	#TODO
	
		###Set initial conditions
	DG_rover.set_initial_values({
		"environment:(temperature)" : 0,
		"thermal interface:(resistance)" : 20000,
		"environment:[dust contamination]" : 0,
		"WEB:(temperature)" : 230,#10+273, #K
		"wheel:(torque)" : 0,
		"wheel:[stall]" : 0,
		"environment:[rock contact]" : 0,
		"wheel:[snag]" : 0,
		"environment:(slope)" : 0,
		"environment:(sinkage)" : 0.02, #2 cm
		"power bus:[transient]" : 0,
		"battery:(power draw)" : 0,
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
		"CONST_wheel_radius" : 0.4, 
		"CONST_wheel_width" : 0.3, 
		"CONST_soil_exponent" : 1.0, 
		"CONST_soil_kc" : 1400, #N/m2
		"CONST_soil_kphi" : 830000, #N/m3 frictional modulus of soil deformation
	})	
	
	
	###Define functions for each node
	#DG_rover.generate_specification_template()
	DG_rover.specify_node("WEB:(temperature)", system_temperature_simulator, {
        "ti_resistance" : "thermal interface:(resistance)",
        "env_temperature": "environment:(temperature)", 
		"battery_state" : "battery:(power draw)",
		"solar_absorptivity" : "CONST_solar_absorptivity", 
		"infrared_emissivity" : "CONST_infrared_emissivity", 
		"surf_area" : "CONST_surface_area"
	})
	DG_rover.specify_node("environment:(temperature)", lunar_temp_simulator)
	DG_rover.specify_node("thermal interface:(resistance)", ti_resistance_simulator, 
	{
        "temp_hist" : "WEB:(temperature)",
        "contaminants" : "environment:[dust contamination]",
		"coeff_thermal_cycling" : "CONST_thermresistance_thermalcycling_dependence",
		"coeff_contamination" : "CONST_thermresistance_contamination_dependence"
	}, get_history=["WEB:(temperature)"])
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
        "wheelsnag_hist" : "wheel:[snag]",
		"lambda_0" : "CONST_wheelstall_rate",
		"snag_dependence" : "CONST_wheelstall_snagdependence",
		"rock_dependence" : "CONST_wheelstall_rockdependence",
		"load_dependence" : "CONST_wheelstall_loaddependence",
	}, get_history=["wheel:(load)","environment:[rock contact]","wheel:[snag]"])
	DG_rover.specify_node("environment:[rock contact]", poisson_process, {
		"lambd" : "CONST_rock_rate"
	})
	DG_rover.specify_node("wheel:[snag]", poisson_process, {
		"lambd" : "CONST_snag_rate"
	})
	DG_rover.specify_node("environment:(slope)", env_slope_sinemodel)
	DG_rover.specify_node("environment:(sinkage)", env_sinkage_sampler, {
		"sinkage_mean" : "CONST_sinkage_mean",
		"sinkage_stddev" : "CONST_sinkage_stddev" 
	})
	DG_rover.specify_node("power bus:[transient]", transient_NHPP, {
        "wheel_load" : "wheel:(load)",
        "system_temp" : "WEB:(temperature)",
        "battery_state" : "battery:(power draw)",
		"max_nominal_temp" : "CONST_max_nominal_temp", 
		"excess_temp_dependence" : "CONST_transient_excesstemp_dependence", 
		"load_dependence" : "CONST_transient_load_dependence"
	})
	DG_rover.specify_node("battery:(power draw)", battery_state_simulator)
	DG_rover.specify_node("power bus:[brownout]", brownout_NHPP, {
        "wheel_load" : "wheel:(load)",
        "system_temp" : "WEB:(temperature)",
        "battery_state" : "battery:(power draw)",
		"max_nominal_temp" : "CONST_max_nominal_temp", 
		"excess_temp_dependence" : "CONST_brownout_excesstemp_dependence", 
		"load_dependence" : "CONST_brownout_load_dependence"
	})
	
	DG_rover.get_specifications()
	
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
	
	if args.run == "extended_problem":
		extended_problem()
		
	if args.run == "restricted_problem":
		simplified_problem()
		
	if args.run == 'debug':
		debug()