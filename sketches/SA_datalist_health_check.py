import os
import sys
import csv

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from itertools import islice

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
	
	###Safely read out all of the samples into matrices
	Ay = []
	By = []
	Cy = []
	
	if doPrint:
		print("Reading the data files...",flush=True)
	
	if do_subset == 0:
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for l,row in enumerate(csvreader):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_A.csv')
				else:
					Ay.append([float(elem) for elem in row])
		
		with open(base_name+'_B.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for l,row in enumerate(csvreader):
				if len(row) != len(var_names)+1:
					if doDiagnostic:
						print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_B.csv')
				else:
					By.append([float(elem) for elem in row])

		for p,name in enumerate(var_names):
			Ciy = []
			with open(base_name+'_C_'+name+'.csv') as csvfile:
				csvreader = csv.reader(csvfile, delimiter=',')
				for l,row in enumerate(csvreader):
					if len(row) != len(var_names)+1:
						if doDiagnostic:
							print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_C_'+name+'.csv')
					else:
						Ciy.append([float(elem) for elem in row])
			Cy.append(Ciy)
	else:
		lim = int(do_subset/2)
		###Optionally, we can analyze less than the full set of provided samples
		with open(base_name+'_A.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
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
				csvreader = csv.reader(csvfile, delimiter=',')
				for l,row in enumerate(islice(csvreader, lim)):
					if len(row) != len(var_names)+1:
						if doDiagnostic:
							print("Warning: dropped line",l+1,"(length "+str(len(row))+')',"from",base_name+'_C_'+name+'.csv')
					else:
						Ciy.append([float(elem) for elem in row])
			Cy.append(Ciy)
			
	return Ay,By,Cy

def temp_cut_rows(Ay,By,lines):
	Ay_cut = []
	By_cut = []
	#Cy_cut = []

	for i,line in enumerate(Ay):
		if i+1 not in lines:
			Ay_cut.append(line)
	
	for i,line in enumerate(By):
		if i+1 not in lines:
			By_cut.append(line)
	
	"""
	for p,name in enumerate(var_names):
		Ciy_cut = []
		for i,line in enumerate(Cy[p]):
			if i+1 not in lines:
				Ciy_cut.append(line)
		Cy_cut.append(Ciy_cut)
	"""
				
	return Ay_cut,By_cut
	
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

if __name__ == '__main__':  
	base_name = "SA_QoI_extra"

	var_names = ["gain_red","gain_gre","gain_blu","rn_red","rn_gre","rn_blu","dc_red","dc_gre","dc_blu","qe_red_prec","qe_gre_prec","qe_blu_prec","vph_red_prec","vph_gre_prec","vph_blu_prec","sl_prec","bg_prec","coll_prec","red_l1_prec","red_l2_prec","red_l3_prec","red_l4_prec","red_l5_prec","red_l6_prec","red_l7_prec","gre_l1_prec","gre_l2_prec","gre_l3_prec","gre_l4_prec","gre_l5_prec","gre_l6_prec","gre_l7_prec","blu_l1_prec","blu_l2_prec","blu_l3_prec","blu_l4_prec","blu_l5_prec","blu_l6_prec","blu_l7_prec","blu_l8_prec","fiber_frd"]
	
	Ay,By,Cy = read_SA_files("SA_QoI_extra", var_names, do_subset=0, doPrint=True)
	Ayo,Byo,Cyo = read_SA_files("SA_QoI", var_names, do_subset=0, doPrint=True)
	Ay.extend(Ayo)
	By.extend(Byo)
	for p,Ciy in enumerate(Cy):
		Cy[p].extend(Cyo[p])
		
	oversave("SA_QoI", var_names, Ay,By,Cy)
			
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
		
		
		