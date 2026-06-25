#This details an obed script that can be called regardless of the problem

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix,precision_score, classification_report, recall_score,f1_score,roc_curve, roc_auc_score
import warnings

sys.path.append('..')
from obed.mcmc import *
from obed.pdf_estimation import *
from inference.goal_based_inference import *
from inference.bn_modeling import *
from approx.kl_divergence import kl_divergence_2bernoullis

"""
This function solves a utility U(d,y,theta) = 1/Var[H(theta|y,d)] using goal-based inference
For each d, we train a GMM on the joint distribution of data y and QoI H
For each data point y in the Monte Carlo, we condition the GMM on the data to estimate the posterior H
And we find the variance of that posterior H, wich is our u that we sum up
"""
	
def U_rd_info_gbi(d, problem, gmm_qyd, n_mc, doPrint=False):   
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	if doPrint:
		print("Generating",n_mc,"joint MC samples of theta and y...",flush=True)
	pthetas = problem.prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	if doPrint:
		print(Y1_list, flush=True)
	
	#We expect we have already trained a joint gmm model p(y,d,q) offline,
	#and will condition it for p(q|y,d)
	if doPrint:
		print("Conditioning GMM in MC loop...",flush=True)
		
	#Pre-calculate an inverse matrix, to speed up the MC loop:
	inv_Sig_dd, logdet_Sig_dd = gbi_precalc_Sigdd(gmm_qyd, p_dim=1)
	
	U_list = []
	for i,y in enumerate(Y1_list): #MC loop		
		vi = np.concatenate((y,d))
	
		#Now, use my posterior predictive to calculate the utility
		#H_var = gbi_var_of_conditional_pp(gmm_qyd, vi, 
		#	inv_Sig_dd_precalc=inv_Sig_dd, logdet_Sig_dd_precalc=logdet_Sig_dd)
		#u = H_var
		u=0
		
		#TODO use conditioned gmm to find posterior of q
		#TODO calculate the criterion
		
		U_list.append(u)
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(u)+'\t', flush=True, end='\r')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U = np.average(U_list)
	return U, U_list
	
def meet_req(Q, r, satisfaction_criterion):
	if satisfaction_criterion == "Q < r":
		R = int(q_train < r)
	elif satisfaction_criterion == "Q <= r":
		R = int(q_train <= r)
	elif satisfaction_criterion == "Q > r":
		R = int(q_train > r)
	elif satisfaction_criterion == "Q >= r":
		R = int(q_train >= r)
	elif satisfaction_criterion == "ra <= Q <= rb":
		ra = r[0]
		rb = r[1]
		R = 1 if ra <= q_train and q_train <= rb else 0
	else:
		#error
		exit
	
def sample_requirements_classifier(problem, sample_file, n_samples, sample_x=False, buffer_rate=1000):
	#First, grab samples
	#(theta, d) -> (y|theta,d, Q|theta)
	#grab them from the bn_sampling method	
	filename = sample_file if sample_file.endswith('.csv') else sample_file+'.csv'
	
	data_buffer = []
	for i in range(n_samples):
		#print("Drawing samples",i,"...",flush=True)
		#Drawing samples theta, d
		theta_sample = problem.prior_rvs(1)
		d_sample = problem.sample_d(1)
		if sample_x:
			x_sample = problem.sample_x(1)
		else:
			x_sample=[]
			
		#Model propagation for Q, y
		q_train = problem.H(theta_sample, x_sample)
		y_train = problem.eta(theta_sample, d_sample)
	
		#Then, obtain the R={0,1} for each
		#R_train = int(problem.check_req(q_train))
		
		#Lastly, save (R, y, d) to file
		save_data = y_train + d_sample + [q_train]
		
		if buffer_rate <= 1: #just stream it into the save file
			with open(filename, 'a+', newline='') as csvfile:
				writer = csv.writer(csvfile)
				writer.writerow(save_data)
		else: #save to buffer, and dump buffer to file at buffer_rate
			data_buffer.append(save_data)
			if (i+1) % buffer_rate == 0:
				print("(dump)",flush=True)
				with open(filename, 'a+', newline='') as csvfile:
					writer = csv.writer(csvfile)
					for data in data_buffer:
						writer.writerow(data)
				data_buffer.clear()
	
def train_requirements_classifier(problem, sample_file):
	#Load (R, y, d) from file
	filename = sample_file if sample_file.endswith('.csv') else sample_file+'.csv'

	###Make sure the files exist
	if not os.path.isfile(filename):
		print("File",filename,"is missing")
		sys.exit()
		
	###Safely read out all of the samples into matrices
	y = []
	d = []
	R = []
	
	print("Reading the data files...",flush=True)
	
	with open(filename) as csvfile:
		csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
		for l,row in enumerate(csvreader):
			if len(row) != problem.dim_y + problem.dim_d + 1:
				if doDiagnostic:
					print("Warning: dropped line",l+1,"(length "+str(len(row))+' expected', str(problem.dim_y + problem.dim_d + 1)+')',"from",filename)
			try:
				ygrab = [float(e) for e in row[:problem.dim_y]]
				dgrab = [float(e) for e in row[problem.dim_y:-1]]
				Qgrab = float(row[-1])
			except ValueError: #recently im seeing some '' values in y? hopefully this avoids that ugliness
				continue
			y.append(ygrab) #should be length dim_y
			d.append(dgrab) #should be length dim_d = row - dim_y - 1
			R.append(int(problem.check_req(Qgrab)))
	n = len(R)
	
	#Train classifier
	yd = [y_i + d_i for y_i,d_i in zip(y,d)]
	yd_names = problem.y_names + problem.d_names
	X = pd.DataFrame(yd, columns=yd_names)
	Y = pd.Series(R)
	
	scaler = StandardScaler()
	X_train = scaler.fit_transform(X)
	#X_test = scaler.transform(X_test)

	classifier = LogisticRegression(max_iter=1000)
	classifier.fit(X_train, Y)
	print("Classifier trained with",n,"samples",flush=True)
	
	#Return classifier
	return classifier, scaler

def validate_requirements_classifier(classifier, scaler, problem, n_val, sample_x=False):
	warnings.filterwarnings("ignore", message="X does not have valid feature names")
	R_true = [0]*n_val
	R_pred = [0]*n_val
	R_predprob = [0]*n_val
	correct = 0
	false_pos = 0
	false_neg = 0
	
	for i in range(n_val):
		#sample a new (theta, d) -> (y|theta,d, Q|theta)
		#evaluate the R_val
		theta_val = problem.prior_rvs(1)
		d_val = problem.sample_d(1)
		x_val = problem.sample_x(1) if sample_x else []
		q_val = problem.H(theta_val, x_val)
		y_val = problem.eta(theta_val, d_val)
		R_val = int(problem.check_req(q_val))
		
		#Evaluate R_class = Class(y_val,d_val)
		yd = np.concatenate((y_val,d_val)).reshape(1,-1)
		xi = scaler.transform(yd)
		R_class = int(classifier.predict(xi)[0])
		R_classprob = classifier.predict_proba(xi)[:, 1]
		
		#spit to screen, with a rolling estimation of the correct, false pos, and false neg
		R_true[i] = R_val
		R_pred[i] = R_class
		R_predprob[i] = R_classprob
		correct += (R_class == R_val)
		false_pos += (R_class > R_val)
		false_neg += (R_class < R_val)
		print("Correct:",correct,"False positive:",false_pos,"False negative:",false_neg, flush=True, end='\r')
		
	print('\n')
	# Output the metrics:
	print(f"Accuracy: {accuracy_score(R_true, R_pred):.2f}")
	print(f"ROC AUC Score: {roc_auc_score(R_true, R_predprob):.2f}")

	# 2. View the Confusion Matrix
	print("\nConfusion Matrix:")
	print(confusion_matrix(R_true, R_pred))

	# 3. View Detailed Precision, Recall, and F1-Scores
	print("\nClassification Report:")
	print(classification_report(R_true, R_pred))
	
#Special approach: solve the requirements-driven problem using a 
#classification-based approach
#First, train a classifier to associate y|d,theta and d with whether 
# Q|theta meets {1} or doesn't meet {0} the requirement r
def U_rd_class(d, problem, n_mc, q_prior, classifier, scaler, doPrint=False):  
	warnings.filterwarnings("ignore", message="X does not have valid feature names") 
	#Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
	if doPrint:
		print("Generating",n_mc,"joint MC samples of theta and y...",flush=True)
	pthetas = problem.prior_rvs(n_mc)
	Y1_list = [problem.eta(theta, d) for theta in pthetas]
	
	#We expect we have already trained a joint gmm model p(y,d,q) offline,
	#and will condition it for p(q|y,d)
	if doPrint:
		print("Running MC loop...",flush=True)
	
	U_list = []
	for i,y in enumerate(Y1_list): #MC loop		
		ydi = np.concatenate((y,d)).reshape(1,-1)
		xi = scaler.transform(ydi)
		
		#use classifier to find the posterior probability of req satisfaction q_posterior
		#R|y,d ~ Ber(q_posterior)
		#Class 0 probability is [0], Class 1 probability is [1]
		q_posterior = classifier.predict_proba(xi)[:, 1]
		#Apply buffers to prevent this from going to 0 or 1, and also strip down to a float
		q_posterior = min(1-1e-11,max(1e-11, q_posterior[0]))
		
		#TODO calculate the criterion
		u_info = kl_divergence_2bernoullis(q_posterior, q_prior)
		u_entropy = (1-q_posterior)*np.log(1-q_posterior) + q_posterior*np.log(q_posterior)
		
		U_list.append([u_info, u_entropy])
		if False:
			print(str(i+1)+"/"+str(n_mc),str(u_info)+'\t'+str(u_entropy), flush=True, end='\r')
		
		if doPrint:
			print(str(i+1)+"/"+str(n_mc),str(q_posterior),str(u_info)+'\t'+str(u_entropy), flush=True, end='\n')
			
	if doPrint:
		print('')
		
	#compute an in-distribution probability
	U_info = np.average([Us[0] for Us in U_list])
	U_entropy = np.average([Us[1] for Us in U_list])
	return U_info, U_entropy, U_list