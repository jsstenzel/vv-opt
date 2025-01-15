#import argparse
import os
import sys
import csv

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

#sys.path.insert(0, "..")
#from problems.problem_definition import *

#Here, I am manually implementing Sobol Saltelli algorithm

"""
Summary of functions:

saltelli_eval
	Take in already-made samples
	construct new A, B, Ci
	evaluate model for all new samples
	append new samples & evaluations onto files
	
saltelli_indices
	Take in filenames
	Calculate and print Si, STi
		
saltelli_eval_sample
	Take in variable names, distributions, bounds/parameters
	Samples from those distributions
	Then calls saltelli_eval

saltelli_eval_sobol
	Take in variable names, distributions, bounds/parameters
	Samples from the SciPy Sobol sequence, starting the sequence in a place that continues any existing Sobol samples
	Performs CDF post-processing on the Sobol sequence to make it match the correct distributions
	Then calls saltelli_eval
	
Example: Saltelli iteration over M
	Take in a list of cumulative sample sizes to use
	Take in a sampler function with 1 argument: sample size
	Function loop:
		Go to each provided sample size
		use saltelli_eval to make up the samples, if necessary
		use saltelli_indices
"""


def saltelli_eval(new_samples, base_name, var_names, model):
	if len(new_samples) % 2 != 0:
		print("Need an even number of new samples to saltelli_eval", flush=True)
		sys.exit()
		
	new_A = []
	new_B = []
	###Split 2M into A and B deterministically
	for m,row in enumerate(new_samples):
		if m % 2 == 0:
			new_A.append(row)
		else:
			new_B.append(row)
			
	###Construct the p Ci matrices
	new_C = []
	for p,name in enumerate(var_names):
		new_Ci = []
		
		#Construct new_Ci by subsituting the pth column of A into B
		for n,b_row in enumerate(new_B):
			a_row = new_A[n]
			
			c_row = b_row[:]
			c_row[p] = a_row[p]
			new_Ci.append(c_row)
		
		#Now save that into a p-long list of C's
		new_C.append(new_Ci)
		
	###Perform the M*(2+p) model evaluations, append onto each row of the new matrices
	for m,a_row in enumerate(new_A):
		yA = model(a_row)
		new_A[m].append(yA)
	
	for m,b_row in enumerate(new_B):
		yB = model(b_row)
		new_B[m].append(yB)
	
	for new_Ci in new_C:
		for m,ci_row in enumerate(new_Ci):
			yCi = model(ci_row)
			new_Ci[m].append(yCi)
		
	###Finally, append these new matrices onto each of the 2+p files
	aw = 'a'# if os.path.exists(base_name+'_A.csv') else 'w+'
	with open(base_name+'_A.csv', aw, newline='') as csvfile:
		writer = csv.writer(csvfile)
		for row in new_A:
			writer.writerow(row)
	
	aw = 'a'# if os.path.exists(base_name+'_B.csv') else 'w+'
	with open(base_name+'_B.csv', aw, newline='') as csvfile:
		writer = csv.writer(csvfile)
		for row in new_B:
			writer.writerow(row)

	for p,new_Ci in enumerate(new_C):
		aw = 'a'# if os.path.exists(base_name+'_C_'+var_names[p]+'.csv') else 'w+'
		with open(base_name+'_C_'+var_names[p]+'.csv', aw, newline='') as csvfile:
			writer = csv.writer(csvfile)
			for row in new_Ci:
				writer.writerow(row)


def saltelli_indices(base_name, var_names, doPrint=True):
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

	###First, read out all of the samples into matrices
	Ay = []
	By = []
	Cy = []
	
	with open(base_name+'_A.csv') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			Ay.append([float(elem) for elem in row])
	
	with open(base_name+'_B.csv') as csvfile:
		csvreader = csv.reader(csvfile, delimiter=',')
		for row in csvreader:
			By.append([float(elem) for elem in row])

	for p,name in enumerate(var_names):
		Ciy = []
		with open(base_name+'_C_'+name+'.csv') as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',')
			for row in csvreader:
				Ciy.append([float(elem) for elem in row])
		Cy.append(Ciy)
		
	###Isolate A,B,Ci and yA,yB,yCi
	A = [Ay_row[:-1] for Ay_row in Ay] #all but last element
	yA = [Ay_row[-1] for Ay_row in Ay] #only last element
	B = [By_row[:-1] for By_row in By]
	yB = [By_row[-1] for By_row in By]
	
	C = []
	yC = []
	for p,Ciy in enumerate(Cy):
		Ci = [Ciy_row[:-1] for Ciy_row in Ciy]
		yCi = [Ciy_row[-1] for Ciy_row in Ciy]
		C.append(Ci)
		yC.append(yCi)

	
	###Provide final calculations
	f02 = np.dot(np.mean(yA),np.mean(yA))  #inner product of p-length and p-length vectors
	S = []
	ST = []
	M = len(yA)
	
	for p,name in enumerate(var_names):
		yAyCi = np.dot(yA, yC[p])  #inner product of n-length and n-length vectors
		yAyA = np.dot(yA, yA)      #inner product of n-length and n-length vectors
		yByCi = np.dot(yB, yC[p])  #inner product of n-length and n-length vectors
		Si = (yAyCi/M - f02)/(yAyA/M - f02)
		STi = 1 - (yByCi/M - f02)/(yAyA/M - f02)
		S.append(Si)
		ST.append(STi)
		
	if doPrint:
		print("-------------------------------------")
		print("Saltelli GSA algorithm, M =",M,"p =",len(var_names))
		print("Total number of model evaluations:",M*(len(var_names)+2))
		print("Parameter",'\t','\t',"S_i",)
		for p,name in enumerate(var_names):
			print(name,'\t','\t',S[p])
		print("Parameter",'\t','\t',"S_Ti")
		for p,name in enumerate(var_names):
			print(name,'\t','\t',ST[p])
		print("-------------------------------------", flush=True)
	
	return S, ST


def saltelli_iterate(sampler, base_name, var_names, model):
	0
	
def clear_saltelli_sample_files(base_name, var_names):
	os.remove(base_name+'_A.csv')
	os.remove(base_name+'_B.csv')
	for p,name in enumerate(var_names):
		os.remove(base_name+'_C_'+var_names[p]+'.csv')


def __ishigami(X):
	a=7
	b=0.1
	x1 = X[0]
	x2 = X[1]
	x3 = X[2]
	return np.sin(x1) + a * np.sin(x2)**2 + b * x3**4 * np.sin(x1)
	
#According to https://uqpyproject.readthedocs.io/en/stable/auto_examples/sensitivity/comparison/ishigami.html#sphx-glr-auto-examples-sensitivity-comparison-ishigami-py
#I should expect the following first-order indices:
#S1 = 0.3139
#S2 = 0.4424
#S3 = 0.0
#and the following total-order indices:
#ST1 = 0.55758886
#ST2 = 0.44241114
#ST3 = 0.24368366

def __print_ishigami_indices(a=7, b=0.1):
	print("Expected first-order indices:")
	VY = a**2/8 + (b*np.pi**4)/5 + b**2 * (np.pi**8)/18 + 0.5 
	V1 = 0.5 * (1 + (b*np.pi**4)/5)**2
	V2 = a**2 / 8
	V3 = 0
	print("S_1 =", V1/VY)
	print("S_2 =", V2/VY)
	print("S_3 =", V3/VY)
	print("Expected total-effect indices:")
	VT1 = 0.5 * (1 + (b*np.pi**4)/5)**2 + (8*b*b*np.pi**8)/225
	VT2 = a**2 / 8
	VT3 = (8*b*b*np.pi**8)/225
	print("S_T1 =", VT1/VY)
	print("S_T2 =", VT2/VY)
	print("S_T3 =", VT3/VY)

if __name__ == '__main__':  
	###########################
	#Full demonstration
	###########################
	
	###Set up the problem
	var_names = ['x1','x2','x3']
	dists = ['norm', 'norm', 'norm']
	bounds = [[-np.pi,np.pi],[-np.pi,np.pi],[-np.pi,np.pi]]
	
	###Generate samples
	samples = []
	for _ in range(10000):
		samples.append(list(scipy.stats.uniform.rvs(size=3, loc=-np.pi, scale=2*np.pi)))
	
	###Test 
	#clear_saltelli_sample_files("ishigami", var_names)
	saltelli_eval(samples, "ishigami", var_names, __ishigami)
	saltelli_indices("ishigami", var_names, doPrint=True)
	
	__print_ishigami_indices()