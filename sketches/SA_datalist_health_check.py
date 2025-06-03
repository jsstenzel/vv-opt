import os
import sys
import csv

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from itertools import islice

#var_names = ["gain_red","gain_gre","gain_blu","rn_red","rn_gre","rn_blu","dc_red","dc_gre","dc_blu","qe_red_prec","qe_gre_prec","qe_blu_prec","vph_red_prec","vph_gre_prec","vph_blu_prec","sl_prec","bg_prec","coll_prec","red_l1_prec","red_l2_prec","red_l3_prec","red_l4_prec","red_l5_prec","red_l6_prec","red_l7_prec","gre_l1_prec","gre_l2_prec","gre_l3_prec","gre_l4_prec","gre_l5_prec","gre_l6_prec","gre_l7_prec","blu_l1_prec","blu_l2_prec","blu_l3_prec","blu_l4_prec","blu_l5_prec","blu_l6_prec","blu_l7_prec","blu_l8_prec","fiber_frd"]
var_names = ["Us","Ud","Qc","I_SMhubt","I_SMhuba","K_yPM","I_xRWA","I_yRWA","I_RWt","I_RWa","I_ISOa","I_ISOt","K_yISO","K_xISO","I_bus","I_propt","I_propa","I_i1","I_i2","I_i3","A_sptop","D_sp","t_sp","I_ss","K_rad1","K_rad2","K_rISO","K_act1","K_act2","I_iso","K_zpet","K_pm1","K_pm3","K_pm4","K_pm5","K_pm6","K_act_pm2","K_act_pm3","K_act_pm4","K_act_pm5","K_act_pm6","K_xpet","c_RWA","c_RWAI","c_SM_act","c_PM","c_PM_act","c_petal","zeta_sunshield","zeta_isolator","zeta_solarpanel","Ru","fc","Tst","Srg","Sst","Tgs","lambda_","Ro","QE","Mgs","fca","Kc","Kcf"]


def read_SA_files(base_name,var_names,do_subset=0,doPrint=False):
	doDiagnostic=doPrint
	###Make sure the files exist
	if not os.path.isfile(base_name+'_A.csv'):
		print("File",base_name+'_A.csv',"is missing")
		sys.exit()
	if not os.path.isfile(base_name+'_B.csv'):
		print("File",base_name+'_B.csv',"is missing")
		sys.exit()
	for p,name in enumerate(var_names):
		if not os.path.isfile(base_name+'_C_'+name+'.csv'):
			print("File",base_name+'_C_'+name+'.csv',"is missing")
			sys.exit()
			
	##check for nulls:
	"""
	if '\0' in open(base_name+'_A.csv').read():
		print("you have null bytes in",base_name+'_A.csv')
		sys.exit()
	if '\0' in open(base_name+'_B.csv').read():
		print("you have null bytes in",base_name+'_B.csv')
		sys.exit()
	for p,name in enumerate(var_names):
		if '\0' in open(base_name+'_C_'+name+'.csv').read():
			print("you have null bytes in",base_name+'_C_'+name+'.csv')
			sys.exit()
	"""
	
	###Safely read out all of the samples into matrices
	Ay = []
	By = []
	Cy = []
	
	if doPrint:
		print("Reading the data files...",flush=True)
	
	if do_subset == 0:
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for l,row in enumerate(csvreader):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_A.csv')
				else:
					Ay.append([float(elem) for elem in row])
		
		with open(base_name+'_B.csv') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for l,row in enumerate(csvreader):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_B.csv')
				else:
					By.append([float(elem) for elem in row])

		for p,name in enumerate(var_names):
			Ciy = []
			print(name)
			with open(base_name+'_C_'+name+'.csv') as csvfile:
				csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
				for l,row in enumerate(csvreader):
					if len(row) != len(var_names)+1:
						if doDiagnostic:
							print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_C_'+name+'.csv')
					else:
						try:
							Ciy.append([float(elem) for elem in row])
						except:
							print(row)
							Ciy.append([0 for elem in row])
			Cy.append(Ciy)
	else:
		lim = int(do_subset/2)
		###Optionally, we can analyze less than the full set of provided samples
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
			for l,row in enumerate(islice(csvreader, lim)):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_A.csv')
				else:
					Ay.append([float(elem) for elem in row])
		
		with open(base_name+'_B.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for l,row in enumerate(islice(csvreader, lim)):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_B.csv')
				else:
					By.append([float(elem) for elem in row])

		for p,name in enumerate(var_names):
			Ciy = []
			with open(base_name+'_C_'+name+'.csv') as csvfile:
				csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
				for l,row in enumerate(islice(csvreader, lim)):
					if len(row) != len(var_names)+1:
						if doDiagnostic:
							print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_C_'+name+'.csv')
					else:
						Ciy.append([float(elem) for elem in row])
			Cy.append(Ciy)
			
	return Ay,By,Cy

def temp_cut_rows(Ay,By,Cy,lines):
	Ay_cut = []
	By_cut = []
	Cy_cut = []

	for i,line in enumerate(Ay):
		if i+1 not in lines:
			Ay_cut.append(line)
	
	for i,line in enumerate(By):
		if i+1 not in lines:
			By_cut.append(line)
	
	for p,name in enumerate(var_names):
		Ciy_cut = []
		for i,line in enumerate(Cy[p]):
			if i+1 not in lines:
				Ciy_cut.append(line)
		Cy_cut.append(Ciy_cut)
				
	return Ay_cut,By_cut,Cy_cut
	
def oversave(base_name, var_names, Ay,By,Cy):
	aw = 'w+'# if os.path.exists(base_name+'_A.csv') else 'w+'
	with open(base_name+'_A.csv', aw, newline='') as csvfile:
		writer = csv.writer(csvfile)
		for row in Ay:
			writer.writerow(row)
	
	aw = 'w+'# if os.path.exists(base_name+'_B.csv') else 'w+'
	with open(base_name+'_B.csv', aw, newline='') as csvfile:
		writer = csv.writer(csvfile)
		for row in By:
			writer.writerow(row)

	for p,Ciy in enumerate(Cy):
		aw = 'w+'# if os.path.exists(base_name+'_C_'+var_names[p]+'.csv') else 'w+'
		with open(base_name+'_C_'+var_names[p]+'.csv', aw, newline='') as csvfile:
			writer = csv.writer(csvfile)
			for row in Ciy:
				writer.writerow(row)

def health_check(base_name, var_names, Ay=[], By=[], Cy=[], do_subset=0):
	if not (Ay and By and Cy):
		print("Health check:",base_name,flush=True)
		Ay,By,Cy = read_SA_files(base_name, var_names, do_subset=do_subset, doPrint=True)
	#else:
		#ignore base name and use the provided lists
			
	###################################################
	#See if A and B match:
	if len(Ay) != len(By):
		print("Big trouble, the file lengths of A and B are different")
		sys.exit()
		
	###################################################
	#See if every row is the expected size:
	for i,line in enumerate(Ay):
		if len(line) != len(var_names)+1:
			print("Length mismatch with A on line",i,"- expected",len(var_names)+1,"got",len(line),flush=True)
	
	for i,line in enumerate(By):
		if len(line) != len(var_names)+1:
			print("Length mismatch with B on line",i,"- expected",len(var_names)+1,"got",len(line),flush=True)
			
	for p,name in enumerate(var_names):
		for i,line in enumerate(Cy[p]):
			if len(line) != len(var_names)+1:
				print("Length mismatch with C_"+name+" on line",i,"- expected",len(var_names)+1,"got",len(line),flush=True)
	
	###################################################
	#See if Ci matches A&B:
	for p,name in enumerate(var_names):
		if len(Ay) != len(Cy[p]):
			print("Trouble, the data lengths of A ("+str(len(Ay))+") and C_"+name+" ("+str(len(Cy[p]))+") are different")

	
	###################################################
	#Try to see if every C is constructed as it ought to be
	
	A = [Ay_row[:-1] for Ay_row in Ay] #all but last element
	B = [By_row[:-1] for By_row in By]	
		
	#Check each parameter to make sure the constructions of the C matrices are right:
	for p,name in enumerate(var_names):
		check_Ci = []
		
		#Construct check_Ci by substituting the pth column of A into B
		for n,b_row in enumerate(B):
			a_row = A[n]
			
			c_row = b_row[:]
			c_row[p] = a_row[p]
			check_Ci.append(c_row)
		
		#Compare check_Ci to Cy:
		first_err = False
		i = 0
		while i < len(check_Ci) and not first_err:
			try:
				if check_Ci[i] != Cy[p][i][:-1]: #comparing two rows:
					print("Mismatch with",name,"on line",str(i+1),flush=True)
					#print(check_Ci[i][0], Cy[p][i][:-1][0])
					first_err = True
			except:
				print("Mismatch with",name,"on last line",flush=True)
				first_err = True
			i += 1
	
	###################################################
	#Lastly, see if there are any duplicates across A+B
	AB_seen = set()
	A_duplicates = []
	B_duplicates = []
	for ii,item in enumerate(A):
		fixed_item = tuple(item)
		if fixed_item in AB_seen:
			A_duplicates.append(ii)
		else:
			AB_seen.add(fixed_item)
	for ii,item in enumerate(B):
		fixed_item = tuple(item)
		if fixed_item in AB_seen:
			B_duplicates.append(ii)
		else:
			AB_seen.add(fixed_item)
	for idx in A_duplicates:
		print("duplicate in A:",idx)
	for idx in B_duplicates:
		print("duplicate in B:",idx, flush=True)

def base_combine(primary_base, var_names, list_bases, save=False):
	Ay,By,Cy = read_SA_files(primary_base, var_names, do_subset=0, doPrint=True)
	
	for base in list_bases:
		Ay2,By2,Cy2 = read_SA_files(base, var_names, do_subset=0, doPrint=True)
		
		#Combine the files
		Ay.extend(Ay2)
		By.extend(By2)
		for p,Ciy in enumerate(Cy2):
			Cy[p].extend(Ciy)

	if save:
		oversave(primary_base, var_names, Ay,By,Cy)
	else:
		health_check("", var_names, Ay,By,Cy)
		
def base_trim(base_name, var_names, do_subset=0):	
	Ay,By,Cy = read_SA_files(base_name, var_names, do_subset=do_subset, doPrint=True)
	oversave(base_name, var_names, Ay,By,Cy)
	
def base_line_cut(base_name, var_names, lines):
	Ay,By,Cy = read_SA_files(base_name, var_names, do_subset=0, doPrint=True)
	temp_cut_rows(Ay,By,Cy,lines)
	oversave(base_name, var_names, Ay,By,Cy)
		
if __name__ == '__main__':

	#base_line_cut("SA_jitter", var_names, [29398,29397])
	health_check("SA_jitter",var_names)
	
	#Check all bases
	#base_names = ["SA_QoI_bank"+str(i) for i in range(1,129)]
	#print(base_names)
	#for base in base_names:
	#	health_check(base, var_names)
		
	#check some trimmed bases
	#health_check("SA_QoI_bank35", var_names, do_subset=246*2)
	#base_trim("SA_QoI_bank35", var_names, do_subset=246*2)
	#base_trim("SA_jitter", var_names, do_subset=29398*2)
	
	#health_check("SA_QoI_bank91", var_names, do_subset=222*2)
	#base_trim("SA_QoI_bank91", var_names, do_subset=222*2)
	
	#See if they all work ok combined
	#base_combine("SA_QoI", var_names, base_names, False)
	
	#combine and save
	#base_combine("SA_QoI", var_names, base_names, True)