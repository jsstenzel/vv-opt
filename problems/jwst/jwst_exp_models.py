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

def apply_error(y, d, scale=0, err=True):
	if err==False:
		return y
		
	if d==0:
		return 0 #the experiment essentially fails, no meaningful data
	elif d==1:
		stddev = scale * 3
	elif d==2:
		stddev = scale
	elif d==3:
		stddev = scale / 3
	
	random = scipy.stats.norm.rvs(scale = stddev)
	meas = y + random
	return meas

def jwst_eta(theta, d, x, prior_mean, err=True):
	#define interest params:
	if type(theta) is not dict:
		print("oh howd that happen now?")
		sys.exit()
	
	#################################################
	###Transmissibility tests: indirectly measure damping coefficient and stiffness of an element
	#Spacecraft Bus - NOT MODELED
	#Spacecraft Bus structure - NOT MODELED
	#RWA Isolator
	y_transmit_rwai = transmissibility_rwaisolator(theta, d["d_transmit_rwai"])
	#CCA Isolator
	y_transmit_ccai = transmissibility_ccaisolator(theta, d["d_transmit_ccai"])
	#Isolator Array (IA)
	y_transmit_ia = transmissibility_isolatorarray(theta, d["d_transmit_ia"])
	#Optics - NOT MODELED
	#Primary Mirrors & Structure
	y_transmit_pmss = transmissibility_primaryassembly(theta, d["d_transmit_pmss"])
	#Secondary Mirrors & Structure
	y_transmit_smss = transmissibility_secondaryassembly(theta, d["d_transmit_smss"])
	#Telescope
	y_transmit_telescope = transmissibility_telescope(theta, d["d_transmit_telescope"])
	
	#################################################
	###Stiffness tests
	
	#RWA Isolator
	y_stiff_xISO = stiffness_measure(theta["K_xISO"], theta["stiffness_rt_factor"], d["d_stiff_rwai"]) #K_xISO
	y_stiff_yISO = stiffness_measure(theta["K_yISO"], theta["stiffness_rt_factor"], d["d_stiff_rwai"]) #K_yISO
	y_stiff_aISO = stiffness_measure((5/3)*theta["K_rISO"], theta["stiffness_rt_factor"], d["d_stiff_rwai"]) #K_aISO
	y_stiff_rISO = stiffness_measure(theta["K_rISO"], theta["stiffness_rt_factor"], d["d_stiff_rwai"]) #K_rISO
	#CCA Isolator - SKIP
	#Primary mirror support structure - NOT MODELED
	#Primary mirrors (deployed)
	y_stiff_pm1 = stiffness_measure(theta["K_pm1"], theta["stiffness_rt_factor"], d["d_stiff_pmss"]) #K_pm1
	y_stiff_yPM = stiffness_measure(theta["K_yPM"], theta["stiffness_rt_factor"], d["d_stiff_pmss"]) #K_yPM
	y_stiff_pm3 = stiffness_measure(theta["K_pm3"], theta["stiffness_rt_factor"], d["d_stiff_pmss"]) #K_pm3
	y_stiff_pm4 = stiffness_measure(theta["K_pm4"], theta["stiffness_rt_factor"], d["d_stiff_pmss"]) #K_pm4
	y_stiff_pm5 = stiffness_measure(theta["K_pm5"], theta["stiffness_rt_factor"], d["d_stiff_pmss"]) #K_pm5
	y_stiff_pm6 = stiffness_measure(theta["K_pm6"], theta["stiffness_rt_factor"], d["d_stiff_pmss"]) #K_pm6
	#Isolator Array (IA) - NOT MODELED
	#Tower Assembly - SKIP
	#Solar array - NOT MODELED
	#Sunshield - NOT MODELED
	#Secondary mirror support structure - NOT MODELED
	#PM Actuator - NOT JWST
	y_stiff_act_pm2 = stiffness_measure(theta["K_act_pm2"], theta["stiffness_rt_factor"], d["d_stiff_pm"]) #K_act_pm2
	y_stiff_act_pm3 = stiffness_measure(theta["K_act_pm3"], theta["stiffness_rt_factor"], d["d_stiff_pm"]) #K_act_pm3
	y_stiff_act_pm4 = stiffness_measure(theta["K_act_pm4"], theta["stiffness_rt_factor"], d["d_stiff_pm"]) #K_act_pm4
	y_stiff_act_pm5 = stiffness_measure(theta["K_act_pm5"], theta["stiffness_rt_factor"], d["d_stiff_pm"]) #K_act_pm5
	y_stiff_act_pm6 = stiffness_measure(theta["K_act_pm6"], theta["stiffness_rt_factor"], d["d_stiff_pm"]) #K_act_pm6
	#Deployable petal
	y_stiff_xpet = stiffness_measure(theta["K_xpet"], theta["stiffness_rt_factor"], d["d_stiff_petal"]) #K_xpet
	y_stiff_zpet = stiffness_measure(theta["K_zpet"], theta["stiffness_rt_factor"], d["d_stiff_petal"]) #K_zpet
	#SM Actuator
	y_stiff_act1 = stiffness_measure(theta["K_act1"], theta["stiffness_rt_factor"], d["d_stiff_sm"]) #K_act1
	y_stiff_act2 = stiffness_measure(theta["K_act2"], theta["stiffness_rt_factor"], d["d_stiff_sm"]) #K_act2
	y_stiff_rad1 = stiffness_measure(theta["K_rad1"], theta["stiffness_rt_factor"], d["d_stiff_sm"]) #K_rad1
	y_stiff_rad2 = stiffness_measure(theta["K_rad2"], theta["stiffness_rt_factor"], d["d_stiff_sm"]) #K_rad2
	
	#################################################
	###Modal surveys: get dominant frequencies and mode shapes, from which we might be able to infer stiffness and transmissibility, but not directly measure
	#Full Observatory - too complicated, SKIP
	
	#Telescope
	y_modal_telescope = transmissibility_telescope(theta, d["d_modal_telescope"])
	#Primary Mirror & Structure
	y_modal_pmss = transmissibility_primaryassembly(theta, d["d_modal_pmss"])
	#Primary Mirror
	#SKIP, not really distinguished for NEXUS
	#Secondary Mirror & Structure
	y_modal_smss = transmissibility_secondaryassembly(theta, d["d_modal_smss"])
	#Optics & Focal Plane -NOT MODELED
	#Spacecraft Bus structure - NOT MODELED
	#Cryo-cooler
	y_modal_ccai = transmissibility_ccaisolator(theta, d["d_modal_ccai"])
	#Sunshield
	y_modal_sunshield = stiffness_measure(theta["zeta_sunshield"], theta["damping_rt_factor"], d["d_modal_sunshield"])
	#Solar Array -- NOT JWST
	y_modal_solar = stiffness_measure(theta["zeta_solarpanel"], theta["damping_rt_factor"], d["d_modal_solar"])
	
	#################################################
	###Micro-Vibe tests: measuring force and moment disturbances generated by these components
	#Reaction wheel disturbance
	#c_RWA reaction wheel damping
	#here, k and m refer to the mass and stiffness of the wheel
	w_wheel = x["Ru"] * (1/60) #rpm -> hz
	y_vibe_rwa = microvibe_RWA(theta["Us"], theta["Ud"], w=w_wheel, c=theta["c_RWA"], k=theta["K_rISO"], m=x["m_RWAchx"], di=d["d_vibe_rwa"])
	
	#Reaction wheel disturbance with isolator
	#c_RWA #reaction wheel damping
	#c_RWAI #reaction wheel isolator damping
	#here, k and m refer to the mass and stiffness of the isolator assembly as a whole
	#use K_rISO and 
	y_vibe_rwai = microvibe_RWA(theta["Us"], theta["Ud"], w=w_wheel, c=theta["c_RWA"]+theta["c_RWAI"], k=theta["K_rISO"], m=x["m_RW"], di=d["d_vibe_rwai"])
	
	#Cryo-cooler disturbance
	y_vibe_cca = microvibe_CCA(C=x["C"], Qc=1, x=x, di=d["d_vibe_cca"]) #TODO make C in theta
	
	#Cryo-cooler disturbance with isolator
	y_vibe_ccai = microvibe_CCA(C=x["C"], Qc=theta["Qc"], x=x, di=d["d_vibe_ccai"])
	
	#################################################
	###Other
	#Cryogenic Modal Survey (characterize shift of stiffness and transmissibility from room temp to cryogenic temp)
	y_cryo_stiffness = simple_measure(theta["stiffness_rt_factor"], d["d_cryo"])
	y_cryo_damping = simple_measure(theta["damping_rt_factor"], d["d_cryo"])
	
	y_full = [
		y_transmit_rwai,
		y_transmit_ccai,
		y_transmit_ia,
		y_transmit_pmss,
		y_transmit_smss,
		y_transmit_telescope,
		y_stiff_xISO,
		y_stiff_yISO,
		y_stiff_aISO,
		y_stiff_rISO,
		y_stiff_pm1,
		y_stiff_yPM,
		y_stiff_pm3,
		y_stiff_pm4,
		y_stiff_pm5,
		y_stiff_pm6,
		y_stiff_act_pm2,
		y_stiff_act_pm3,
		y_stiff_act_pm4,
		y_stiff_act_pm5,
		y_stiff_act_pm6,
		y_stiff_xpet,
		y_stiff_zpet,
		y_stiff_act1,
		y_stiff_act2,
		y_stiff_rad1,
		y_stiff_rad2,
		y_modal_telescope,
		y_modal_pmss,
		y_modal_smss,
		y_modal_ccai,
		y_modal_sunshield,
		y_modal_solar,
		y_vibe_rwa,
		y_vibe_rwai,
		y_vibe_cca,
		y_vibe_ccai,
		y_cryo_stiffness,
		y_cryo_damping
	]
	return y_full

###############################################################################################
###############################################################################################
#Definitions & full info ######################################################################
###############################################################################################
###############################################################################################


#################################################
###Transmissibility tests: indirectly measure damping coefficient and stiffness of an element
#Spacecraft Bus
"""
This would include the reaction wheels and everything else above the IA. However for NEXUS model, nothing else in the bus is really modeled.
NOT MODELED which is unfortunate but I don't really have a better idea other than making the model quite a bit more detailed..."""

#Spacecraft Bus structure
"""NOT MODELED really at all for same reasons as above"""

#RWA Isolator
"""related to nicelas2 rows 141 to 218
K_xISO, K_yISO, K_aISO, K_rISO
ADDED c_RWA, but thats just for the wheels, not the isolator?
There are these other points as well - not RWA attach points, but with the same stiffness properties. 
I think ill assume those are isolator assembly points, ADDED c_RWAI
"""
def transmissibility_rwaisolator(theta, di):
	#K_xISO, K_yISO, K_aISO, K_rISO, c_RWAI
	K_xISO = theta["K_xISO"] * theta["stiffness_rt_factor"]
	K_yISO = theta["K_yISO"] * theta["stiffness_rt_factor"]
	c_RWAI = theta["c_RWAI"] * theta["damping_rt_factor"]
	
	y = c_RWAI / math.sqrt(K_xISO**2 + K_xISO**2 + K_yISO**2)
	
	meas = apply_error(y, di)
	return meas

#CCA Isolator
"""Node 207 is defined to be where the vibration enters in. This is in the definition of ig near line 2027
(Note that this is the only isolation between the instrument and the cryocooler for JWST, they are separated by the isolator array, but not here)
(Could model the instrument point in nicelas2 with stiffness K_cryo and damping c_cryo, but don't want to hack the NEXUS model representation)
there is no CCA isolator in the NEXUS design
However, the "cryocooler attenuation factor" Qc effectively models this I think
"""
def transmissibility_ccaisolator(theta, di):
	Qc = theta["Qc"] * theta["damping_rt_factor"]
	
	y = Qc
	
	meas = apply_error(y, di)
	return meas

#Isolator Array (IA)
"""On JWST, there is a single structure that bridges the bus to the rest of the spacecraft
On NEXUS, it is a little less clear. I think 64,27,28,30 are the points that connect telescope structure to bus
However, they are not modeled in stiffness matrix
(Could model the 4 connection points in nicelas2 with stiffness K_IA and damping c_IA, treat that as IA equivalent structure, but don't want to hack the NEXUS model representation)
The effectively equivalent structure is the isolators on the support structure
"""
def transmissibility_isolatorarray(theta, di):
	zeta_isolator = theta["zeta_isolator"] * theta["damping_rt_factor"]

	y = zeta_isolator
	
	meas = apply_error(y, di)
	return meas

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
def transmissibility_primaryassembly(theta, di):
	K_pm1 = theta["K_pm1"] * theta["stiffness_rt_factor"]
	K_yPM = theta["K_yPM"] * theta["stiffness_rt_factor"]
	K_pm3 = theta["K_pm3"] * theta["stiffness_rt_factor"]
	K_act2 = theta["K_act2"] * theta["stiffness_rt_factor"]
	K_act_pm2 = theta["K_act_pm2"] * theta["stiffness_rt_factor"]
	K_act_pm3 = theta["K_act_pm3"] * theta["stiffness_rt_factor"]
	K_xpet = theta["K_xpet"] * theta["stiffness_rt_factor"]
	K_zpet = theta["K_zpet"] * theta["stiffness_rt_factor"]
	c_PM = theta["c_PM"] * theta["damping_rt_factor"]
	c_PM_act = theta["c_PM_act"] * theta["damping_rt_factor"]

	#K_pm1, K_yPM, K_pm3, K_pm4, K_pm5, K_pm6, c_PM
	y_PM = c_PM / math.sqrt(K_pm1**2 + K_yPM**2 + K_pm3**2)
	#K_act2, K_act_pm2, K_act_pm3, K_act_pm4, K_act_pm5, K_act_pm6, c_PM_act
	y_PM_act = c_PM_act / math.sqrt(K_act2**2 + K_act_pm2**2 + K_act_pm3**2)
	#K_xpet, K_zpet, c_petal
	y_xpet = c_PM_act / math.sqrt(K_xpet**2 + K_zpet**2 + K_xpet**2)
	y = y_PM * y_PM_act * y_xpet
	
	meas = apply_error(y, di)
	return meas

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
def transmissibility_secondaryassembly(theta, di):
	#K_act1, K_act2, K_rad1, K_rad2, c_SM_act
	K_act1 = theta["K_act1"] * theta["stiffness_rt_factor"]
	K_act2 = theta["K_act2"] * theta["stiffness_rt_factor"]
	c_SM_act = theta["c_SM_act"] * theta["damping_rt_factor"]
	
	y = c_SM_act / math.sqrt(K_act1**2 + K_act2**2 + K_act2**2)
	
	meas = apply_error(y, di)
	return meas

#Telescope
"""Includes everything downstream from the 4 bus connection points
Primary mirrors and structure + secondary mirrors and structure
"""	
def transmissibility_telescope(theta, di):
	K_pm1 = theta["K_pm1"] * theta["stiffness_rt_factor"]
	K_yPM = theta["K_yPM"] * theta["stiffness_rt_factor"]
	K_pm3 = theta["K_pm3"] * theta["stiffness_rt_factor"]
	K_act2 = theta["K_act2"] * theta["stiffness_rt_factor"]
	K_act_pm2 = theta["K_act_pm2"] * theta["stiffness_rt_factor"]
	K_act_pm3 = theta["K_act_pm3"] * theta["stiffness_rt_factor"]
	K_xpet = theta["K_xpet"] * theta["stiffness_rt_factor"]
	K_zpet = theta["K_zpet"] * theta["stiffness_rt_factor"]
	K_act1 = theta["K_act1"] * theta["stiffness_rt_factor"]
	c_SM_act = theta["c_SM_act"] * theta["damping_rt_factor"]
	c_PM = theta["c_PM"] * theta["damping_rt_factor"]
	c_PM_act = theta["c_PM_act"] * theta["damping_rt_factor"]
	zeta_isolator = theta["zeta_isolator"] * theta["damping_rt_factor"]
	
	#K_pm1, K_yPM, K_pm3, K_pm4, K_pm5, K_pm6, c_PM
	y_PM = c_PM / math.sqrt(K_pm1**2 + K_yPM**2 + K_pm3**2)
	#K_act2, K_act_pm2, K_act_pm3, K_act_pm4, K_act_pm5, K_act_pm6, c_PM_act
	y_PM_act = c_PM_act / math.sqrt(K_act2**2 + K_act_pm2**2 + K_act_pm3**2)
	#K_xpet, K_zpet, c_petal
	y_xpet = c_PM_act / math.sqrt(K_xpet**2 + K_zpet**2 + K_xpet**2)
	#K_act1, K_act2, K_rad1, K_rad2, c_SM_act
	y_SM_act = c_SM_act / math.sqrt(K_act1**2 + K_act2**2 + K_act2**2)
	#zeta_isolator
	y = y_PM * y_PM_act * y_xpet * y_SM_act * zeta_isolator
	
	meas = apply_error(y, di)
	return meas


#################################################
###Stiffness tests
def stiffness_measure(K, stiffness_factor, di):
	y = K * stiffness_factor
	
	meas = apply_error(y, di)
	return meas

#RWA Isolator
"""Directly measure K_xISO, K_yISO, K_aISO, K_rISO"""
#K_xISO
#K_yISO
#K_aISO
#K_rISO

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
#K_pm1
#K_yPM
#K_pm3
#K_pm4
#K_pm5
#K_pm6

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

#Sunshield
"""NOT MODELED - at least, not directly as stiffness
There is a modal damping for the solar array that is modeled, kind of an analog
"""

#Secondary mirror support structure
"""The secondary mirror support structure doesnt have stiffness defined - theyre assumed to be perfectly rigid
NOT MODELED"""

#PM Actuator
#This wasnt one of the tests identified for JWST, but its relevant for this modeling approach...?
#K_act2
#K_act_pm2
#K_act_pm3
#K_act_pm4
#K_act_pm5
#K_act_pm6

#Deployable petal
#K_xpet
#K_zpet

#SM Actuator
#K_act1
#K_act2
#K_rad1
#K_rad2


#################################################
###Modal surveys: get dominant frequencies and mode shapes, from which we might be able to infer stiffness and transmissibility, but not directly measure
###Note that these experiments are outside the jitter scope strictly, and they're done anyway for launch survivability testing
#Full Observatory

#Telescope
#just use transmissibility_telescope with higher uncertainty

#Primary Mirror & Structure
#just use transmissibility_primaryassembly with higher uncertainty

#Primary Mirror
#SKIP, not really distinguished for NEXUS

#Secondary Mirror & Structure
#just use transmissibility_secondaryassembly with higher uncertainty

#Optics & Focal Plane
#NOT MODELED

#Spacecraft Bus structure
#NOT MODELED

#Cryo-cooler
#just use transmissibility_ccaisolator with higher uncertainty

#Sunshield
#simple_measure zeta_sunshield

#Solar Array -- this one is not historical
#simple_measure zeta_solarpanel


#################################################
###Micro-Vibe tests: measuring force and moment disturbances generated by these components
#Reaction wheel disturbance
#model this really simply, where both US and Ud are creating a force magnitude
#Us [gcm]
#Ud [gcm^2]
#c_RWA #reaction wheel damping
#here, k and m refer to the mass and stiffness of the wheel
def microvibe_RWA(Us, Ud, w, c, k, m, di):
	reaction_wheel_radius = 0.393 / 2
	static_force = Us * w * w #force = static imbalance * angularmomentum^2
	dynamic_torque = Ud * w * w #torque = dynamic imbalance * angularmomentum^2
	dynamic_force = dynamic_torque / reaction_wheel_radius
	F = static_force + dynamic_force #magnitude of external force on system
	
	#I care about the amplitude of the forced spring-mass-damping system:
	#https://courses.lumenlearning.com/suny-osuniversityphysics/chapter/15-6-forced-oscillations
	#Amax = F0/sqrt(m(w^2 - w0^2) + c^2 w^2)
	w0 = math.sqrt(k/m)
	y = F / math.sqrt(m*(w**2 - w0**2)**2 + c**2 * w**2)
	
	meas = apply_error(y, di)
	return meas

#Reaction wheel disturbance with isolator
#now apply damping to that foce magnitude
#I think same as above, but I want to see the disturbance of the larger system
#so now, ill add the two damping values:
#c_RWA #reaction wheel damping
#c_RWAI #reaction wheel isolator damping
#here, k and m refer to the mass and stiffness of the isolator assembly as a whole
#use K_rISO and 

#Cryo-cooler disturbance
def microvibe_CCA(C, Qc, x, di):
	#x:
	#n=3
	#h = [1.0,2.0,3.0,4.0,5.0,6.0]
	#C = [
	#		[42,0.95,4.1,2.75,0.9,1.2],
	#		[0.2,0.09,0.25,1.0,5.0,0.4]
	#	]
	#fc = 30 #ish
	n = x["n"]
	h = x["h"]
	#C = x["C"] #maybe i should make this a theta? otherwise nothing is tested here
	fc = x["fc"]
	wc=2*math.pi*fc # drive frequency

	#do smarter indexing to do the right thing here
	zij = np.zeros((2,n))
	kij = np.zeros((2,n))
	for ind1 in [0,1]:
		for ind2 in range(n):
			zij[ind1][ind2]=1/(2*wc*h[ind2])
			kij[ind1][ind2]=2*zij[ind1][ind2]*wc*h[ind2]*Qc*C[ind1][ind2]
	print(zij)
	print(kij)

	#fix this
	Adc=[ #NxN = 12x12
			[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[-wc**2, -2*zij[0][0]*wc, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, -2**2*wc**2, -2*zij[0][1]*2*wc, 0, 0, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, -3**2*wc**2, -2*zij[0][2]*3*wc, 0, 0, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, -1**2*wc**2, -2*zij[1][0]*1*wc, 0, 0, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, -2**2*wc**2, -2*zij[1][1]*2*wc, 0, 0],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
			[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3**2*wc**2, -2*zij[1][2]*3*wc],
	]
	Bdc=[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1] #NxP = 12x1
	Cdc=[ #QxN = 3x12
			[0, 0, 0, 0, 0, 0, 0, kij[1][0], 0, kij[1][1], 0, kij[1][2]],
			[0, 0, 0, 0, 0, 0, 0, kij[1][0], 0, kij[1][1], 0, kij[1][2]],
			[0, kij[0][0], 0, kij[0][1], 0, kij[0][2], 0, 0, 0, 0, 0, 0]
		]
	Ddc=[0, 0, 0] #QxP = 3x1
	Adc = np.array(Adc)
	Bdc = np.array(Bdc).reshape(12,1)
	Cdc = np.array(Cdc)
	Ddc = np.array(Ddc).reshape(3,1)
	
	#do state space solver and see what values we get
	cryo = scipy.signal.StateSpace(Adc, Bdc, Cdc, Ddc)
	#cryo.to_discrete(t) #use the same timestep as the simulation
	ts = np.arange(0, 30, 0.1) #time vector
	u = np.array([[1] for e in ts]) #input vector, no idea what i would use here
	x0 = [0]*12 #initial conditions, no idea what i would use here
	#print(u)
	#print(ts)
	#print(x0)
	t, y, x = scipy.signal.lsim(cryo, U=u, T=ts, X0=x0)
	#print(t) # = len(ts)
	#print(y) #Q = 3
	#print(x) #N = 12
	#plt.plot(y)
	#plt.show()
	
	#test this to extract a kind of disturbance measure from it
	y_peak = min(yi[2] for yi in y)
	
	#TODO add error with di
	meas = apply_error(y_peak, di)
	return meas

#Cryo-cooler disturbance with isolator
#call microvibe_cryocooler with a nonzero Qc

#################################################
###Other
#Cryogenic Modal Survey (characterize shift of stiffness and transmissibility from room temp to cryogenic temp)
#simple measure, include stiffness_rt_factor and damping_rt_factor in theta, put into every experiment
def simple_measure(theta, di):
	y = theta
	
	meas = apply_error(y, di)
	return meas