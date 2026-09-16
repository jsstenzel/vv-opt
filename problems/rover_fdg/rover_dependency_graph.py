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
	DG_rover = SystemDependencyGraph(timestep=0.1)
	
	###MISSION: identify quantities of interest
	DG_rover.add("mission:(downtime)")
	DG_rover.add("mission:(productive traversal)",'<-',["control:(operational mode)", "vehicle:(average drive speed)"]) #drive at speed if in drive mode
	
	###FAULTS: think about the failure modes that could affect QoI
	
	#whenever there is unproductive traversal for whatever reason, it prevents that much productive traversal. i.e. if i would move 10km, but 4km of that is detour, then i have made 6km productive traversal
	DG_rover.add("mission:(unproductive traversal remaining)",'<->',"mission:(productive traversal)")
	DG_rover.add("mission:(unproductive traversal remaining)",'<-',["control:[avoid obstacle]", "software:[nav error]"])
	
	#transient and brownout can trigger a software/hardware reset
	DG_rover.add("power bus:[transient]","->","mission:(downtime)")
	DG_rover.add("power bus:[brownout]","->","mission:(downtime)")
	
	#a stall or snag leads to a brief downtime
	DG_rover.add("actuator:[stall]","->","mission:(downtime)")
	DG_rover.add("wheel:[snag]","->","mission:(downtime)")
	
	#loss of an actuator causes a decrease in average speed
	DG_rover.add("actuator:[loss]","->","control:(average drive speed)")
	
	#accumulation of motor drive electronics stress causes a decrease in average speed
	DG_rover.add("motor drive electronics:[fault]",'->',"control:(average drive speed)")
	DG_rover.add("motor drive electronics:[fault]",'<-',"motor drive electronics:(stress)")
	
	#deterioration of the thermal interface & thermal resistance increasing leads to uncontrolled temperature
	#this leads to more cooldown stops, greater fault probabilities due to temperature
	DG_rover.add("thermal interface:(resistance)",'->',"WEB:(temperature)")
	DG_rover.add("WEB:(temperature)","->","control:(operational mode)")
	
	#could add wheel degradation here
	
	###FDIR: think about how observing, recovering, and mitigating those faults can introduce additional effects
	
	#detection&mitigation of actuator loss involves changing traversal pattern, partially resolving loss of speed but increasing wear on remaining actuators
	DG_rover.add("fdir:[actuator loss detected]",'<-',"actuator:[loss]") #detection
	DG_rover.add("fdir:[actuator loss mitigated]",'<-',"fdir:[actuator loss detected]") #detection
	DG_rover.add("fdir:[actuator loss mitigated]",'->',"actuator:[loss]") #exacerbation
	DG_rover.add("fdir:[actuator loss mitigated]",'->',"control:(average drive speed)") #mitigation
	
	#a detected motor drive electionics fault requires drive abort and stuck recovery
	DG_rover.add("fdir:[mde fault detected]",'<-',"motor drive electronics:[fault]")
	DG_rover.add("fdir:[mde fault mitigated]",'<-',"fdir:[mde fault detected]")
	DG_rover.add("fdir:[mde fault mitigated]",'->',"mission:(downtime)")
	
	###SYSTEM & CONOPS
	
	#pose is determined nominally by # productive traversal, but detours change that
	DG_rover.add("vehicle:(pose)","<-",["mission:(productive traversal)", "mission:(unproductive traversal remaining)"])
	
	#sleep mode at night, safe mode if we're incurring downtime
	DG_rover.add("control:(operational mode)","<-",["environment:(UTC time)", "mission:(downtime)"])
	
	DG_rover.add("control:(average drive speed)") #with each actuator down, we have to go slower
	
	###ENVIRONMENTAL EFFECTS
	
	#greater dust contamination during the day
	DG_rover.add("environment:[dust contamination]","<-","environment:(temperature)")
	
	#model temperature due to time of lunar day
	DG_rover.add("environment:(temperature)","<-",["environment:(UTC time)"])
	
	#slip is a property of soil at a given location, partially determined by the slope and sinkage at that location
	DG_rover.add("environment:(slip ratio)","<-",["vehicle:(pose)", "environment:(slope)", "environment:(sinkage)"])
	DG_rover.add("environment:(slope)","<-",["vehicle:(pose)"])
	DG_rover.add("environment:(sinkage)","<-",["vehicle:(pose)", "environment:(slope)"])
	
	#obstacles may include large rocks, areas of unfavorable terrain
	DG_rover.add("control:[avoid obstacle]","<-",["environment:[obstacle]"])
	
	###THERMAL
	
	DG_rover.add("radiator:(efficiency)","<-",["environment:[dust contamination]"])
	DG_rover.add("thermal interface:(resistance)","<-",["WEB:(temperature)", "environment:[dust contamination]"])
	DG_rover.add("WEB:(temperature)","<-",["thermal interface:(resistance)", "environment:(temperature)", "battery:(power draw)", "radiator:(efficiency)"])
	
	###MOBILITY
	
	#stall occurs when torque hits a limit, and that limit is lower if mde stress is high. MDE stress accumulates with torque
	DG_rover.add("actuator:[stall]","<-",["wheel:(torque)", "motor drive electronics:(stress)"]) #discrete event, when torque hits a limit
	DG_rover.add("motor drive electronics:(stress)","<-",["wheel:(torque)"])
	
	#will lose actuator according to torque exposure, with higher probability as stalls and dust acucmulate
	DG_rover.add("actuator:[loss]","<-",["wheel:(torque)", "actuator:[stall]", "environment:[dust contamination]"])
	
	#many instantaneous determinants of torque
	DG_rover.add("wheel:(torque)","<-",["environment:(slope)", "environment:(slip ratio)", "environment:(sinkage)", "environment:[rock contact]", "wheel:[snag]", "wheel:(degradation)", "control:(average drive speed)", "actuator:[loss]"]) #continuous variable
	
	#wheel degrades with all distance ever traversed, and rocks and snags
	DG_rover.add("wheel:(degradation)","<-",["mission:(distance traversed)", "environment:[rock contact]", "wheel:[snag]"])
	DG_rover.add("mission:(distance traversed)","<-",["mission:(productive traversal)", "mission:(unproductive traversal remaining)"]) #simple addition
	DG_rover.add("wheel:[snag]","<-",["environment:[rock contact]"])
	
	#driving over high-slip terrain kicks up dust, requires us to slow down
	DG_rover.add("environment:(slip ratio)",'->',"environment:[dust contamination]")
	DG_rover.add("environment:(slip ratio)",'->',"control:(average drive speed)")
	
	###POWER
	
	DG_rover.add("power bus:[transient]","<-",["wheel:(torque)", "WEB:(temperature)"]) #discrete event
	DG_rover.add("power bus:[brownout]","<-",["WEB:(temperature)", "battery:(power draw)"]) #discrete event
	DG_rover.add("battery:(power draw)","<-","control:(operational mode)")
	
	###SOFTWARE
	
	DG_rover.add("software:(estimated pose)","<-","vehicle:(pose)")
	DG_rover.add("software:[nav error]",'<-',["software:(estimated pose)", "vehicle:(pose)"])
	
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
	DG_rover.printout()
	
	###Set categories and display
	DG_rover_categorization = {
		"control:(fault recovery)" : "CONTROL",
		"mission:(downtime)" : "MISSION",
		"mission:(productive traversal)" : "MISSION",
		"mission:(distance traversed)" : "MISSION",
		"mission:(unproductive traversal remaining)" : "FAULT",
		"control:(operational mode)" : "CONTROL",
		"control:(average drive speed)" : "CONTROL",
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
		"fdir:[actuator loss detected]" : "OBSERVATION",
		"fdir:[actuator loss mitigated]" : "CONTROL",
		"fdir:[mde fault detected]" : "OBSERVATION",
		"fdir:[mde fault mitigated]" : "CONTROL",
		"vehicle:(pose)" : "MISSION",
		"environment:(UTC time)" : "ENVIRONMENT",
		"environment:[obstacle]" : "FAULT",
		"environment:(temperature)" : "ENVIRONMENT",
		"environment:(slip ratio)" : "ENVIRONMENT",
		"environment:(slope)" : "ENVIRONMENT",
		"environment:(sinkage)" : "ENVIRONMENT",
		"radiator:(efficiency)" : "SYSTEM",
		"thermal interface:(resistance)" : "SYSTEM",
		"wheel:(degradation)" : "SYSTEM",
		"software:[nav error]" : "FAULT",
		"software:(estimated pose)" : "OBSERVATION",
	}
	DG_rover.categorize_nodes(DG_rover_categorization)
	
	category_colors = {
		"FAULT":"red",
		"QOI":"black",
		"ENVIRONMENT":"green",
		"SYSTEM":"orange",
		"OBSERVATION":"blue",
		"CONTROL":"purple",
	}
	DG_rover.set_category_colors(category_colors)
	
	DG_rover.display_interactive()

	
def restricted_problem():
	###Define graph structure
	DG_rover = SystemDependencyGraph(timestep=0.1)
	
	###MISSION: identify quantities of interest
	DG_rover.add("mission:(downtime)")
	#productive traversal is a subset of total traversal
	DG_rover.add("mission:(distance traversed)","->",["mission:(productive traversal)"]) #simple addition
	DG_rover.add("mission:(distance traversed)",'<-',["control:(operational mode)", "vehicle:(average drive speed)"]) #drive at speed if in drive mode


	###FAULTS and first-order effects
	DG_rover.add("vehicle:(average drive speed)",'<-',"control:(average drive speed)") #drive at speed if in drive mode
	
	#a stall or snag leads to a brief downtime
	DG_rover.add("actuator:[stall]","->","mission:(downtime)")
	
	#accumulation of motor drive electronics stress causes an involuntary decrease in average speed
	DG_rover.add("motor drive electronics:[fault]",'->',"vehicle:(average drive speed)")
	DG_rover.add("wheel:(torque)",'->',"motor drive electronics:[fault]")
	

	###FDIR: think about how observing, recovering, and mitigating those faults can introduce additional effects
	
	#a detected motor drive electionics fault requires drive abort and stuck recovery
	#this resolves the slowdown
	DG_rover.add("fdir:[mde fault detected]",'<-',"motor drive electronics:[fault]")
	DG_rover.add("fdir:[mde fault mitigated]",'<-',"fdir:[mde fault detected]")
	DG_rover.add("fdir:[mde fault mitigated]",'->',"mission:(downtime)")
	DG_rover.add("fdir:[mde fault mitigated]",'->',"vehicle:(average drive speed)")
	"""
	###SYSTEM & CONOPS
	#pose is determined nominally by # productive traversal, but detours change that
	DG_rover.add("vehicle:(pose)","<-",["mission:(productive traversal)"])
	
	#sleep mode at night, safe mode if we're incurring downtime
	DG_rover.add("control:(operational mode)","<-",["environment:(UTC time)", "mission:(downtime)"])
	
	###ENVIRONMENTAL EFFECTS
	
	#slip is a property of soil at a given location, partially determined by the slope and sinkage at that location
	DG_rover.add("environment:(slip ratio)","<-","vehicle:(pose)")
	DG_rover.add("environment:(slope)","<-","vehicle:(pose)")
	DG_rover.add("environment:(sinkage)","<-","vehicle:(pose)")
	#DG_rover.add("environment:(sinkage)","<-","environment:(slope)") #no, its actually that the slope increases slip which increases sinkage
	#DG_rover.add("environment:(slip ratio)","<-","environment:(slope)") #i dont want to tangle with this yet; terramechanics

	#Assertion: when the slip increases, the wheel digs into the regolith, leading to increasing sinkage
	DG_rover.add("environment:(slip ratio)","->","environment:(sinkage)")

	###MOBILITY
	
	#stall occurs when torque hits a limit, and that limit is lower if mde stress is high. MDE stress accumulates with torque
	DG_rover.add("actuator:[stall]","<-",["wheel:(torque)", "motor drive electronics:[fault]"]) #discrete event, when torque hits a limit
	DG_rover.add("motor drive electronics:[fault]","<-",["wheel:(torque)"])
	
	#many instantaneous determinants of torque
	DG_rover.add("wheel:(torque)","<-",["environment:(slope)", "environment:(slip ratio)", "environment:(sinkage)", "environment:[rock contact]", "wheel:(degradation)", "vehicle:(average drive speed)"]) #continuous variable
	#wheel:(deformation)
	
	#wheel degrades with all distance ever traversed, and rocks and snags
	DG_rover.add("wheel:(degradation)","<-",["mission:(distance traversed)", "environment:[rock contact]"])
	DG_rover.add("environment:[rock contact]","<-","mission:(distance traversed)")	
	
	#driving over high-slip terrain kicks up dust, requires us to slow down
	DG_rover.add("environment:(slip ratio)",'->',"control:(average drive speed)")
	"""

	"""
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
	})"""
	DG_rover.printout()
	
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
	
	"""
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
	"""
	
		###Set categories and display
	DG_rover_categorization = {
		"control:(fault recovery)" : "CONTROL",
		"mission:(downtime)" : "QOI",
		"mission:(productive traversal)" : "QOI",
		"mission:(distance traversed)" : "SYSTEM",
		"mission:(unproductive traversal remaining)" : "FAULT",
		"control:(operational mode)" : "CONTROL",
		"control:(average drive speed)" : "CONTROL",
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
		"fdir:[actuator loss detected]" : "OBSERVATION",
		"fdir:[actuator loss mitigated]" : "CONTROL",
		"fdir:[mde fault detected]" : "OBSERVATION",
		"fdir:[mde fault mitigated]" : "CONTROL",
		"vehicle:(pose)" : "SYSTEM",
		"vehicle:(average drive speed)" : "SYSTEM",
		"environment:(UTC time)" : "ENVIRONMENT",
		"environment:[obstacle]" : "FAULT",
		"environment:(temperature)" : "ENVIRONMENT",
		"environment:(slip ratio)" : "ENVIRONMENT",
		"environment:(slope)" : "ENVIRONMENT",
		"environment:(sinkage)" : "ENVIRONMENT",
		"radiator:(efficiency)" : "SYSTEM",
		"thermal interface:(resistance)" : "SYSTEM",
		"wheel:(degradation)" : "FAULT",
		"software:[nav error]" : "FAULT",
		"software:(estimated pose)" : "OBSERVATION",
	}
	DG_rover.categorize_nodes(DG_rover_categorization)
	
	category_colors = {
		"FAULT":"red",
		"QOI":"black",
		"ENVIRONMENT":"green",
		"SYSTEM":"orange",
		"OBSERVATION":"blue",
		"CONTROL":"purple",
	}
	DG_rover.set_category_colors(category_colors)
	
	DG_rover.display_interactive()
	
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
		restricted_problem()
		
	if args.run == 'debug':
		debug()