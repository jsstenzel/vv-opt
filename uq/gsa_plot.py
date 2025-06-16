import os
import sys
import csv

import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

def plot_gsa_full(varnames, S1, ST, S1_conf=[], ST_conf=[], title="", coplot=False, screening=0, xspin=True):
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
		if xspin:
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
		if xspin:
			plt.xticks(rotation=90)
		plt.tight_layout()
		#if logplot:
		#	plt.yscale('log')	
		plt.title(title)
			
		#add whiskers for conf interval	
		if ST_conf:
			plt.errorbar(varnames, ST, yerr=ST_conf, fmt='|', color='r')
		plt.grid(axis = 'y')
		
		if screening:
			plt.axhline(y=screening, color='lightblue', linestyle='-', linewidth=2)
		
		plt.show()
	
	else:
		###Make two adjacent subplots sharing an axis
		fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
		plt.xlabel('Parameters', fontweight ='bold', fontsize = 10)
		if xspin:
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
		["Us",-0.0209,0.0027,0.0321,0.0059],
		["Ud",0.1590,0.0025,0.2016,0.0054],
		["Qc",-0.0205,0.0027,0.0318,0.0059],
		["I_SMhubt",-0.0215,0.0027,0.0315,0.0059],
		["I_SMhuba",-0.0215,0.0027,0.0315,0.0059],
		["K_yPM",-0.0215,0.0027,0.0320,0.0059],
		["I_xRWA",-0.0215,0.0027,0.0315,0.0059],
		["I_yRWA",-0.0170,0.0027,0.0385,0.0059],
		["I_RWt",-0.0187,0.0027,0.0358,0.0059],
		["I_RWa",-0.0215,0.0027,0.0315,0.0059],
		["I_ISOa",-0.0216,0.0027,0.0316,0.0059],
		["I_ISOt",-0.0193,0.0027,0.0325,0.0059],
		["K_yISO",-0.0216,0.0027,0.0317,0.0059],
		["K_xISO",-0.0215,0.0027,0.0315,0.0059],
		["I_bus",-0.0167,0.0027,0.0374,0.0059],
		["I_propt",-0.0217,0.0027,0.0317,0.0059],
		["I_propa",-0.0215,0.0027,0.0315,0.0059],
		["I_i1",-0.0216,0.0027,0.0316,0.0059],
		["I_i2",-0.0205,0.0027,0.0321,0.0059],
		["I_i3",-0.0215,0.0027,0.0315,0.0059],
		["A_sptop",-0.0216,0.0027,0.0318,0.0059],
		["D_sp",0.0438,0.0026,0.1019,0.0057],
		["t_sp",0.0150,0.0027,0.0713,0.0058],
		["I_ss",-0.0215,0.0027,0.0315,0.0059],
		["K_rad1",-0.0215,0.0027,0.0315,0.0059],
		["K_rad2",-0.0215,0.0027,0.0315,0.0059],
		["K_rISO",0.5188,0.0020,0.6030,0.0044],
		["K_act1",-0.0215,0.0027,0.0315,0.0059],
		["K_act2",-0.0215,0.0027,0.0315,0.0059],
		["I_iso",-0.0215,0.0027,0.0315,0.0059],
		["K_zpet",0.0905,0.0026,0.1189,0.0056],
		["K_pm1",-0.0215,0.0027,0.0315,0.0059],
		["K_pm3",-0.0215,0.0027,0.0315,0.0059],
		["K_pm4",-0.0215,0.0027,0.0315,0.0059],
		["K_pm5",-0.0215,0.0027,0.0315,0.0059],
		["K_pm6",-0.0215,0.0027,0.0315,0.0059],
		["K_act_pm2",-0.0217,0.0027,0.0317,0.0059],
		["K_act_pm3",-0.0215,0.0027,0.0315,0.0059],
		["K_act_pm4",-0.0215,0.0027,0.0315,0.0059],
		["K_act_pm5",-0.0215,0.0027,0.0315,0.0059],
		["K_act_pm6",-0.0216,0.0027,0.0315,0.0059],
		["K_xpet",-0.0215,0.0027,0.0315,0.0059],
		["c_RWA",-0.0215,0.0027,0.0315,0.0059],
		["c_RWAI",-0.0215,0.0027,0.0315,0.0059],
		["c_SM_act",-0.0215,0.0027,0.0315,0.0059],
		["c_PM",-0.0216,0.0027,0.0316,0.0059],
		["c_PM_act",-0.0215,0.0027,0.0315,0.0059],
		["c_petal",0.0245,0.0027,0.0824,0.0058],
		["zeta_sunshield",-0.0215,0.0027,0.0315,0.0059],
		["zeta_isolator",-0.0215,0.0027,0.0315,0.0059],
		["zeta_solarpanel",-0.0215,0.0027,0.0316,0.0059],
		["Ru",-0.0215,0.0027,0.0315,0.0059],
		["fc",-0.0215,0.0027,0.0315,0.0059],
		["Tst",-0.0215,0.0027,0.0315,0.0059],
		["Srg",-0.0215,0.0027,0.0315,0.0059],
		["Sst",-0.0215,0.0027,0.0315,0.0059],
		["Tgs",-0.0215,0.0027,0.0315,0.0059],
		["lambda_",-0.0215,0.0027,0.0315,0.0059],
		["Ro",-0.0215,0.0027,0.0315,0.0059],
		["QE",-0.0215,0.0027,0.0315,0.0059],
		["Mgs",-0.0215,0.0027,0.0315,0.0059],
		["fca",-0.0215,0.0027,0.0315,0.0059],
		["Kc",-0.0215,0.0027,0.0315,0.0059],
		["Kcf",-0.0215,0.0027,0.0315,0.0059]
	]
	
	varnames = [row[0] for row in data]
	S1 = [max(0,row[1]) for row in data]
	S1_conf = [row[2] for row in data]
	ST = [max(0,row[3]) for row in data]
	ST_conf = [row[4] for row in data]
	plot_gsa_full(varnames, S1, ST, S1_conf=S1_conf, ST_conf=ST_conf, coplot=False, screening=0.05)

	