import os
import sys
import csv

import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

def plot_gsa_full(varnames, S1, ST, S1_conf=[], ST_conf=[], title="", coplot=False, screening=0):
	# set width of bar
	barWidth = 0.3
	barColor = "lightgrey"
	edgeColor = "grey"

	if coplot == False:
		##########################################
		# Make the S1 plot
		plt.bar(varnames, S1, color=barColor, width = barWidth, edgecolor=edgeColor,)
		
		# Adding Xticks
		plt.xlabel('Parameters', fontweight ='bold', fontsize = 10)
		plt.ylabel("Sobol' main-effect index", fontweight ='bold', fontsize = 10)
		#plt.xticks([r + barWidth for r in range(len(varnames))], varnames)
		plt.xticks(rotation=90)
		plt.tight_layout()
		#if logplot:
		#	plt.yscale('log')	
		plt.title(title)
			
		#add whiskers for conf interval	
		if S1_conf:
			plt.errorbar(varnames, S1, yerr=S1_conf, fmt='|', color='r')
		plt.grid(axis = 'y')
		
		plt.show()
		
		##########################################
		# Make the ST plot
		plt.bar(varnames, ST, color=barColor, width = barWidth, edgecolor=edgeColor,)
		
		# Adding Xticks
		plt.xlabel('Parameters', fontweight ='bold', fontsize = 10)
		plt.ylabel("Sobol' total-order index", fontweight ='bold', fontsize = 10)
		#plt.xticks([r + barWidth for r in range(len(varnames))], varnames)
		plt.xticks(rotation=90)
		plt.tight_layout()
		#if logplot:
		#	plt.yscale('log')	
		plt.title(title)
			
		#add whiskers for conf interval	
		if ST_conf:
			plt.errorbar(varnames, ST, yerr=ST_conf, fmt='|', color='r')
		plt.grid(axis = 'y')
		
		plt.show()
	
	else:
		###Make two adjacent subplots sharing an axis
		fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
		plt.xlabel('Parameters', fontweight ='bold', fontsize = 10)
		plt.xticks(rotation=90)
		plt.tight_layout()
		plt.title(title)
		
		###S1
		ax1.bar(varnames, S1, color=barColor, width = barWidth, edgecolor=edgeColor,)
		ax1.set_ylabel("Sobol' main-effect index", fontweight ='bold', fontsize = 10)
		if S1_conf:
			ax1.errorbar(varnames, S1, yerr=S1_conf, fmt='|', color='r')
		ax1.grid(axis = 'y')
		
		###ST
		ax2.bar(varnames, ST, color=barColor, width = barWidth, edgecolor=edgeColor,)
		ax2.set_ylabel("Sobol' total-order index", fontweight ='bold', fontsize = 10)
		if ST_conf:
			ax2.errorbar(varnames, ST, yerr=ST_conf, fmt='|', color='r')
		ax2.grid(axis = 'y')
		
		if screening:
			ax2.axhline(y=screening, color='lightblue', linestyle='-', linewidth=2)
		
		plt.subplots_adjust(hspace=0)
		plt.show()
	
if __name__ == '__main__':  
	###########################
	#Full demonstration
	###########################
	
	data = [
		#Var name,    S       S_conf          ST     ST_conf
		#["gain_red",    -0.0464, 0.0070,  0.0969,  0.0140],
		#["gain_gre",    -0.0465, 0.0070,  0.0969,  0.0140],
		#["gain_blu",    -0.0467, 0.0070,  0.0972,  0.0140],
		["rn_red",      -0.0456, 0.0070,  0.0987,  0.0140],
		["rn_gre",      -0.0448, 0.0070,  0.0968,  0.0140],
		["rn_blu",      -0.0468, 0.0070,  0.0991,  0.0140],
		["dc_red",      0.1894,  0.0063,  0.3126,  0.0128],
		["dc_gre",      0.0989,  0.0065,  0.2825,  0.0132],
		["dc_blu",      0.1414,  0.0064,  0.2646,  0.0130],
		["qe_red_prec", -0.0353, 0.0070,  0.1029,  0.0140],
		["qe_gre_prec", -0.0488, 0.0070,  0.1178,  0.0140],
		["qe_blu_prec", 0.0059,  0.0068,  0.1240,  0.0138],
		["vph_red_prec",0.0528,  0.0066,  0.1956,  0.0134],
		["vph_gre_prec",0.1401,  0.0064,  0.2362,  0.0131],
		["vph_blu_prec",0.0339,  0.0068,  0.1521,  0.0137],
		["sl_prec",     -0.0461, 0.0070,  0.0966,  0.0140],
		["bg_prec",     -0.0463, 0.0070,  0.0969,  0.0140],
		["coll_prec",   -0.0467, 0.0070,  0.0999,  0.0140],
		["red_l1_prec", -0.0473, 0.0070,  0.0980,  0.0140],
		["red_l2_prec", -0.0464, 0.0070,  0.0970,  0.0140],
		["red_l3_prec", -0.0468, 0.0070,  0.0974,  0.0140],
		["red_l4_prec", -0.0468, 0.0070,  0.0974,  0.0140],
		["red_l5_prec", -0.0465, 0.0070,  0.0970,  0.0140],
		["red_l6_prec", -0.0466, 0.0070,  0.0971,  0.0140],
		["red_l7_prec", -0.0461, 0.0070,  0.0968,  0.0140],
		["gre_l1_prec", -0.0464, 0.0070,  0.0969,  0.0140],
		["gre_l2_prec", -0.0463, 0.0070,  0.0968,  0.0140],
		["gre_l3_prec", -0.0470, 0.0070,  0.0976,  0.0140],
		["gre_l4_prec", -0.0466, 0.0070,  0.0972,  0.0140],
		["gre_l5_prec", -0.0464, 0.0070,  0.0969,  0.0140],
		["gre_l6_prec", -0.0455, 0.0070,  0.0967,  0.0140],
		["gre_l7_prec", -0.0470, 0.0070,  0.0975,  0.0140],
		["blu_l1_prec", -0.0465, 0.0070,  0.0971,  0.0140],
		["blu_l2_prec", -0.0463, 0.0070,  0.0968,  0.0140],
		["blu_l3_prec", -0.0461, 0.0070,  0.0966,  0.0140],
		["blu_l4_prec", -0.0464, 0.0070,  0.0969,  0.0140],
		["blu_l5_prec", -0.0468, 0.0070,  0.0973,  0.0140],
		["blu_l6_prec", -0.0463, 0.0070,  0.0969,  0.0140],
		["blu_l7_prec", -0.0464, 0.0070,  0.0969,  0.0140],
		["blu_l8_prec", -0.0463, 0.0070,  0.0968,  0.0140],
		["fiber_frd",   0.0015,  0.0069 , 0.1221,  0.0139]
	]
	
	varnames = [row[0] for row in data]
	S1 = [max(0,row[1]) for row in data]
	S1_conf = [row[2] for row in data]
	ST = [max(0,row[3]) for row in data]
	ST_conf = [row[4] for row in data]
	plot_gsa_full(varnames, S1, ST, S1_conf=S1_conf, ST_conf=ST_conf, coplot=True, screening=0.05)

	