#import argparse
import os
import sys
import shutil
import csv
import fileinput


def cost_factor(d):
	if d==0:
		return 0
	elif d==1:
		return 10000 #$10K
	elif d==2:
		return 100000 #$100K
	elif d==3:
		return 10000000 #$10M

def jwst_cost(d, x):
		cost = 0
		
		component = 1
		assembly = 3
		subsystem = 6
		system = 10
		
		#transmissibility tests
		cost += cost_factor(d["d_transmit_rwai"]) * assembly
		cost += cost_factor(d["d_transmit_ccai"]) * assembly
		cost += cost_factor(d["d_transmit_ia"]) * assembly
		cost += cost_factor(d["d_transmit_pmss"]) * subsystem
		cost += cost_factor(d["d_transmit_smss"]) * subsystem
		cost += cost_factor(d["d_transmit_telescope"]) * system
		
		#stiffness tests
		cost += cost_factor(d["d_stiff_rwai"]) * component
		cost += cost_factor(d["d_stiff_pmss"]) * component
		cost += cost_factor(d["d_stiff_pm"]) * component
		cost += cost_factor(d["d_stiff_petal"]) * component
		cost += cost_factor(d["d_stiff_sm"]) * component
		
		#modal surveys
		cost += cost_factor(d["d_modal_telescope"]) * system
		cost += cost_factor(d["d_modal_pmss"]) * subsystem
		cost += cost_factor(d["d_modal_smss"]) * subsystem
		cost += cost_factor(d["d_modal_ccai"]) * assembly
		cost += cost_factor(d["d_modal_sunshield"]) * assembly
		cost += cost_factor(d["d_modal_solar"]) * assembly
		
		#microvibe tests
		cost += cost_factor(d["d_vibe_rwa"]) * component
		cost += cost_factor(d["d_vibe_rwai"]) * assembly
		cost += cost_factor(d["d_vibe_cca"]) * component
		cost += cost_factor(d["d_vibe_ccai"]) * assembly
		
		#cryogenic materials test
		cost += cost_factor(d["d_cryo"]) * subsystem
		
		return cost