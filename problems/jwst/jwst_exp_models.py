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
	TODO SKIP which is unfortunate but I don't really have a better idea other than making the model quite a bit more detailed..."""
	
	#Spacecraft Bus structure
	"""Not really modeled at all for same reasons as above, SKIP"""
	
	#RWA Isolator
	"""related to nicelas2 rows 141 to 218
	K_xISO, K_yISO, K_aISO, K_rISO
	TODO Need to add 1 damping coefficient c_RWA"""
	
	#CCA Isolator
	"""Node 207 is defined to be where the vibration enters in. This is in the definition of ig near line 2027
	(Note that this is the only isolation between the instrument and the cryocooler; for JWST, they are separated by the isolator array, but not here)
	TODO Need to model the instrument point in nicelas2 with stiffness K_cryo and damping c_cryo
	"""
	
	#Isolator Array (IA)
	"""On JWST, there is a single structure that bridges the bus to the rest of the spacecraft
	On NEXUS, it is a little less clear. I think 64,27,28,30 are the points that connect telescope structure to bus
	However, they are not modeled in stiffness matrix
	TODO model the 4 connection points in nicelas2 with stiffness K_IA and damping c_IA, treat that as IA equivalent structure"""
	
	#Optics
	"""The optics are not modeled in NEXUS, SKIP"""
	
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
	TODO Need to implement c_PM damping coefficient for the actuators 228-401
	"""
	
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
	TODO Need to implement c_SM damping coefficient for the actuators (21-93)
	"""
	
	#Telescope
	"""Includes everything downstream from the 4 bus connection points
	Primary mirrors and structure + secondary mirrors and structure
	"""
	
	#################################################
	###Stiffness tests
	#RWA Isolator
	"""Directly measure K_xISO, K_yISO, K_aISO, K_rISO"""
	
	#CCA Isolator
	"""There does not seem to be any such structure... but the instrument point is defined as the location where cryocooler vibration enters in
	TODO Need to model the instrument point in nicelas2 with stiffness K_cryo"""
	
	#Primary mirror support structure
	"""The primary mirror support structure doesnt have spring-stiffness defined - theyre assumed to be perfectly rigid
	TODO SKIP unless i need to model more things..."""
	
	#Primary mirrors (deployed)
	"""The mirrors (29,115,122,123,124,125,130,140,141,142,143,144,145,160,162,163,164,165,170,180,181,182,183,184,185,190,191,192,193,199,200,201) have stiffness defined:
		0.10000E+07 -> TODO K_pm1
		K_yPM
		0.58400E+06 -> TODO K_pm3
		0.59820E+02 -> TODO K_pm4
		0.49000E+02 -> TODO K_pm5
		0.33250E+02 -> TODO K_pm6
	"""
	
	#Isolator Array (IA)
	"""There is not currently an isolator array.
	TODO model the 4 connection points in nicelas2 with stiffness K_IA, treat that as IA equivalent structure"""
	
	#Tower Assembly
	"""No such structure, theres no tower, just a direct join to bus. SKIP"""
	
	#Solar array
	"""TODO"""
	
	#Sunshield
	"""TODO"""
	
	#Secondary mirror support structure
	"""The secondary mirror support structure doesnt have stiffness defined - theyre assumed to be perfectly rigid
	TODO SKIP unless i need to model more things..."""
	
	
	
	#################################################
	###Modal surveys: get dominant frequencies and mode shapes, from which we might be able to infer stiffness and transmissibility, but not directly measure
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