#import argparse
import os
import sys
import shutil
import csv
import fileinput

import matplotlib.pyplot as plt
from astropy.io import fits
import scipy

import numpy as np
import itertools
import multiprocessing as mp
import math
from copy import deepcopy
import scipy.optimize as optimization

"""
Full matrix experiment model
"""
def jwst_eta(theta, d, x, prior_mean, err=True):
	#define interest params:
	if type(theta) is not dict:
		print("oh howd that happen now?")
		sys.exit()
	#define design variables, gain rn dc:
	t_gain = d["t_gain"]
	
	#################################################
	###Transmissibility tests: indirectly measure damping coefficient and stiffness of an element
	#Spacecraft Bus
	"""
	This would include the reaction wheels and everything else above the IA. However for NEXUS model, nothing else in the bus is really modeled.
	NOT MODELED which is unfortunate but I don't really have a better idea other than making the model quite a bit more detailed..."""
	
	#Spacecraft Bus structure
	"""Not really modeled at all for same reasons as above, SKIP"""
	
	#RWA Isolator
	"""related to nicelas2 rows 141 to 218
	K_xISO, K_yISO, K_aISO, K_rISO
	ADDED c_RWA, but thats just for the wheels, not the isolator?
	There are these other points as well - not RWA attach points, but with the same stiffness properties. 
	I think ill assume those are isolator assembly points, ADDED c_RWAI
	"""
	K_xISO, K_yISO, K_aISO, K_rISO, c_RWAI
	
	#CCA Isolator
	"""Node 207 is defined to be where the vibration enters in. This is in the definition of ig near line 2027
	(Note that this is the only isolation between the instrument and the cryocooler; for JWST, they are separated by the isolator array, but not here)
	(Could model the instrument point in nicelas2 with stiffness K_cryo and damping c_cryo, but don't want to hack the NEXUS model representation)
	there is no CCA isolator in the NEXUS design
	However, the "cryocooler attenuation factor" Qc effectively models this I think
	"""
	Qc
	
	#Isolator Array (IA)
	"""On JWST, there is a single structure that bridges the bus to the rest of the spacecraft
	On NEXUS, it is a little less clear. I think 64,27,28,30 are the points that connect telescope structure to bus
	However, they are not modeled in stiffness matrix
	(Could model the 4 connection points in nicelas2 with stiffness K_IA and damping c_IA, treat that as IA equivalent structure, but don't want to hack the NEXUS model representation)
	The effectively equivalent structure is the isolators on the support structure
	"""
	zeta_isolator
	
	#Optics
	"""NOT MODELED"""
	
	#Primary Mirrors & Structure
	"""3 primary mirrors, telescope support structure
	The mirrors (29,115,122,123,124,125,130,140,141,142,143,144,145,160,162,163,164,165,170,180,181,182,183,184,185,190,191,192,193,199,200,201) have stiffness defined:
		0.10000E+07 -> TODO K_pm1
		K_yPM
		0.58400E+06 -> TODO K_pm3
		0.59820E+02 -> TODO K_pm4
		0.49000E+02 -> TODO K_pm5
		0.33250E+02 -> TODO K_pm6
	The actuators (126,127,128,146,147,148,166,167,168) have stiffnesses defined: 
		K_act2
		0.29100E+07 -> TODO K_act_pm2
		0.10000E+07 -> TODO K_act_pm3
		0.33250E+02 -> TODO K_act_pm4
		0.49000E+02 -> TODO K_act_pm5
		0.12012E+03 -> TODO K_act_pm6
	The primary mirror support structure doesnt have spring-stiffness defined - theyre assumed to be perfectly rigid
	TODO Need to implement c_PM_act damping coefficient for the actuators 228-401
	TODO need to implement c_PM
	TODO need to implement K_xpet = 0.18000E+10, c_petal
	"""
	K_pm1, K_yPM, K_pm3, K_pm4, K_pm5, K_pm6, c_PM
	K_act2, K_act_pm2, K_act_pm3, K_act_pm4, K_act_pm5, K_act_pm6, c_PM_act
	K_xpet, c_petal
	
	#Secondary Mirrors & Structure
	"""secondary mirror, spider
	The mirror (116,117,118,119,120,121,202) dont have spring-stiffness defined.
	The actuators (14,15,16,17,18,19,20,21,23,24,25,26) have stiffnesses defined: 
		K_act1, 
		K_act2, 
		K_act2, 
		K_rad1, 
		K_rad2
		K_rad2
	The secondary mirror support structure doesnt have stiffness defined - theyre assumed to be perfectly rigid
	TODO Need to implement c_SM_act damping coefficient for the actuators (21-93)
	"""
	K_act1, K_act2, K_rad1, K_rad2, c_SM_act
	
	#Telescope
	"""Includes everything downstream from the 4 bus connection points
	Primary mirrors and structure + secondary mirrors and structure
	"""
	K_pm1, K_yPM, K_pm3, K_pm4, K_pm5, K_pm6, c_PM
	K_act2, K_act_pm2, K_act_pm3, K_act_pm4, K_act_pm5, K_act_pm6, c_PM_act
	K_xpet, c_petal
	K_act1, K_act2, K_rad1, K_rad2, c_SM_act
	zeta_isolator
	
	
	#################################################
	###Stiffness tests
	#RWA Isolator
	"""Directly measure K_xISO, K_yISO, K_aISO, K_rISO"""
	K_xISO, K_yISO, K_aISO, K_rISO
	
	#CCA Isolator
	"""There does not seem to be any such structure... but the instrument point is defined as the location where cryocooler vibration enters in
	SKIP this is not modeled with stiffness, but instead as a kind of directly-attenuated system?"""
	
	#Primary mirror support structure
	"""The primary mirror support structure doesnt have spring-stiffness defined - theyre assumed to be perfectly rigid
	NOT MODELED"""
	
	#Primary mirrors (deployed)
	"""The mirrors (29,115,122,123,124,125,130,140,141,142,143,144,145,160,162,163,164,165,170,180,181,182,183,184,185,190,191,192,193,199,200,201) have stiffness defined:
		0.10000E+07 -> TODO K_pm1
		K_yPM
		0.58400E+06 -> TODO K_pm3
		0.59820E+02 -> TODO K_pm4
		0.49000E+02 -> TODO K_pm5
		0.33250E+02 -> TODO K_pm6
	"""
	K_pm1, K_yPM, K_pm3, K_pm4, K_pm5, K_pm6
	
	#Isolator Array (IA)
	"""There is not currently an isolator array.
	(could model the 4 connection points in nicelas2 with stiffness K_IA, treat that as IA equivalent structure)
	The isolators on the primary structure are an analogous structure, but they're NOT MODELED
	There is a modal damping for the isolators that is modeled, kind of an analog
	But since stiffness is NOT MODELED, and I already have a zeta_isolator experiment, SKIP
	"""
	
	
	#Tower Assembly
	"""No such structure, theres no tower, just a direct join to bus. SKIP"""
	
	#Solar array
	"""NOT MODELED - at least, not directly as stiffness
	There is a modal damping for the solar array that is modeled, kind of an analog
	"""
	zeta_solarpanel
	
	#Sunshield
	"""NOT MODELED - at least, not directly as stiffness
	There is a modal damping for the solar array that is modeled, kind of an analog
	"""
	zeta_isolator
	
	#Secondary mirror support structure
	"""The secondary mirror support structure doesnt have stiffness defined - theyre assumed to be perfectly rigid
	NOT MODELED"""
	zeta_solarpanel
	
	
	
	#################################################
	###Modal surveys: get dominant frequencies and mode shapes, from which we might be able to infer stiffness and transmissibility, but not directly measure
	###Note that these experiments are outside the jitter scope strictly, and they're done anyway for launch survivability testing
	#Full Observatory
	#Telescope
	#Primary Mirror & Structure
	#Primary Mirror
	#Secondary Mirror & Structure
	#Optics & Focal Plane
	#Spacecraft Bus structure
	#Cryo-cooler
	#Sunshield
	
	#################################################
	###Micro-Vibe tests: measuring force and moment disturbances generated by these components
	#Reaction wheel disturbance
	#Reaction wheel disturbance with isolator
	#Cryo-cooler disturbance
	#Cryo-cooler disturbance with isolator
	
	#################################################
	###Other
	#Cryogenic Modal Survey (characterize shift of stiffness and transmissibility from room temp to cryogenic temp)



def transmissibility_test():
	0
	#I think I had a different good idea for this test: Measure damping ratio ζ of a number of components, apply noise. ζ is proportional to damping coefficient and invprop to square root of stiffness. In fact ζ = c / 2sqrt(k*m) for simple one-spring-mass-damp system. But how do I actually translate a larger structure into that??
	
	#I think I need to test this by performing the K & M calculation in the NEXUS code,
	#but on a smaller subsystem? That lets me arrive at a set of modal K and M matrices,
	#from which I can extract y measurements??
	#TODO pick up on these ideas


def stiffness_test():
	return stiffness + err