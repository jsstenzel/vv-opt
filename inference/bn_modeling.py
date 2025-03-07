#Building on L&W 2014 and goal_based_inference.py
#Here, I am creating methods for easy and efficient offline training of a GMM(y, d, Q)
#And the methods to sample iteratively, save and load samples, and show convergence of training
#The idea is that, at the end, I can create a GMM definition that itself can be saved and loaded (a list of means, vars, and coefficients)
#Compare to the general approach in saltelli_gsa

import os
import sys
import math
import csv
import numpy as np
from sklearn.mixture import GaussianMixture
import pickle

sys.path.append('..')
from inference.goal_based_inference import *
from uq.plotmatrix import *

#Draw N samples and save in a consistent 1-line format to savefile
#Do it iteratively, to support parallelization and clustering
def bn_sampling(problem, savefile, N, buffer_rate=1, doPrint=False):
	filename = savefile if savefile.endswith('.csv') else savefile+'.csv'
	
	data_buffer = []
	for i in range(N):
		if doPrint:
			print("Drawing Q & y samples",i,"...",flush=True)
		#Drawing samples theta, d
		theta_sample = problem.prior_rvs(1)
		d_sample = problem.sample_d(1)
			
		#Model propagation for Q, y
		qoi_train = problem.H(theta_sample)
		y_train = problem.eta(theta_sample, d_sample)
		
		#Append the new BN sample to file
		save_data = y_train + d_sample + [qoi_train]
		
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
				

#Read and interpret the savefile
def bn_load_samples(problem, savefile, doPrint=False, do_subset=0, doDiagnostic=False):
	filename = savefile if savefile.endswith('.csv') else savefile+'.csv'

	###Make sure the files exist
	if not os.path.isfile(filename):
		print("File",filename,"is missing")
		sys.exit()
		
	###Safely read out all of the samples into matrices
	y = []
	d = []
	Q = []
	
	if doPrint:
		print("Reading the data files...",flush=True)
	
	with open(filename) as csvfile:
		csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
		for l,row in enumerate(csvreader):
			if len(row) != problem.dim_y + problem.dim_d + 1:
				if doDiagnostic:
					print("Warning: dropped line",l+1,"(length "+str(len(row))+' expected', str(problem.dim_y + problem.dim_d + 1)+')',"from",filename)
			else:
				ygrab = [float(e) for e in row[:problem.dim_y]]
				dgrab = [float(e) for e in row[problem.dim_y:-1]]
				Qgrab = float(row[-1])
				y.append(ygrab) #should be length dim_y
				d.append(dgrab) #should be length dim_d = row - dim_y - 1
				Q.append(Qgrab)
	
	#zip it together at the end, in case i ever need them separate for something later
	yd = [y_i + d_i for y_i,d_i in zip(y,d)]
	
	if do_subset:
		Q = Q[:do_subset]
		yd = yd[:do_subset]

	return Q, yd

#Call bn_load_samples and gbi_train_model
def bn_train_from_file(problem, savefile, do_subset=0, doPrint=False):
	###Load file
	qoi_train, y_d_train = bn_load_samples(problem, savefile, doPrint, do_subset)
	#y_d_train = [[yd[0], yd[1]] for yd in y_d_train] #stupid cut for speed
	
	###look at the covariance matrix and eigenvalues of the data, see if we have bad correlation
	if False:
		data_together = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
		names_together = ["Q"] + problem.y_names + problem.d_names
		cov = covmatrix_heatmap(data_together, names_together, rescale=False)
		eigenval, _ = np.linalg.eig(cov)
		print("eigenvalues of the Qyd covariance matrix:",eigenval)
	
	###Train model	
	gmm = gbi_train_model(qoi_train, y_d_train, verbose=2, ncomp=0, careful=True)
	
	###Print and return
	if doPrint:
		print("Trained GMM with",len(qoi_train),"samples from",savefile)
		print(gmm,len(gmm.means_[0]),"dimensions", flush=True)
		
	return gmm
	
#This calculates MSE between the truth and model
#Note that of course, conditioning the GMM gets you the posterior predictive of Q
#So what I actually want to compare here is how well Var[GMM(y,d)] compares to Var[Q]
#NOTE that this doesn't compare Var[GMM(y,d)] to Var[p(Q|y)], to do that I need to do inference, such as with MCMC
#This is a little busted; we're comparing Q data prior to Q model posterior which is apples and pears
#TODO someday implement MCMC in here to compare data MCMC posterior to model GBI posterior
def bn_measure_model_mse(problem, gmm, N, doPrint=True):
	N_val = N if N>0 else 1000
	if doPrint:
		print("Drawing",N_val,"validation samples...",flush=True)
	
	###Draw N_val totally new validation samples
	theta_samples = problem.prior_rvs(N_val)
	d_samples = problem.sample_d(N_val)
	
	###################################
	if doPrint:
		print("Evaluating truth values...",flush=True)
	###Evaluate the validation set for "truth" case
	#for each theta,d get a bunch of subsamples of Q, and calculate the variance of that
	subsample_N = N_val #i suppose?
	Qvar_true = []
	for i in range(N_val):
		qoi_train = [problem.H(theta_samples[i]) for i in range(subsample_N)]
		#y_train = [problem.eta(theta_samples[i], d_samples[i]) for _ in range(subsample_N)]
		Qvar_true.append(np.var(qoi_train, ddof=0))  #not ddof=1, since I want to compare with the conditioned GMM, which uses the "actual" variance!
	
	###################################
	if doPrint:
		print("Evaluating model values...",flush=True)
	###Evaluate the validation set with the GMM
	y_samples = [problem.eta(theta_samples[j], d_samples[j]) for j in range(N_val)]
	#this introduces uncertainty that is not present in the truth case. To calculate real MSE, i would need to actually calculate posterior
		
	v = [y+d for y,d in zip(y_samples, d_samples)] #turn into array
	Qvar_gmm = [gbi_var_of_conditional_pp(gmm, vi, verbose=2) for vi in v] #turn into array
	
	SE = [(truth - model)**2 for truth,model in zip(Qvar_true, Qvar_gmm)]
	MSE = np.mean(SE)
	if doPrint:
		print("MSE for",N_val,"validation samples:",MSE, flush=True)
		
	return MSE

#Save the GMM to file for easy grabbing later
def bn_save_gmm(gmm, gmm_file):
	filename = gmm_file if gmm_file.endswith('.pkl') else gmm_file+'.pkl'
	
	with open(filename, 'wb') as file:
		pickle.dump(gmm, file)
	
	print("GMM saved to",filename,flush=True)

#Load GMM from file
def bn_load_gmm(gmm_file):
	filename = gmm_file if gmm_file.endswith('.pkl') else gmm_file+'.pkl'
	
	with open(filename, 'rb') as file:
		data_gmm = pickle.load(file)

	return data_gmm

def bn_measure_stability_convergence(problem, big_savefile, N_val, doPrint=True):
	N_val = N_val if N_val>0 else 5		
	N_list = [1000,4000,10000,40000,100000,400000]
	
	###Draw N_val totally new validation samples
	#TODO i probably should not be reusing validation samples.
	if doPrint:
		print("Drawing",N_val,"validation samples...",flush=True)
	theta_samples = problem.prior_rvs(N_val)
	d_samples = problem.sample_d(N_val)
	y_samples = [problem.eta(theta_samples[j], d_samples[j]) for j in range(N_val)]
	
	###Use the savedata to train a series of GMMs, save them if they don't already exist; or, load them if they do
	list_gmm = []
	for N in N_list:
		gmmfile_n = "BN_model_" + str(N) + '.pkl'
		if os.path.exists(gmmfile_n):
			###Load it
			if doPrint:
				print("Loading",gmmfile_n,"...",flush=True)
			gmm_n = bn_load_gmm(gmmfile_n)
			list_gmm.append(gmm_n)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",gmmfile_n,"...",flush=True)
			gmm_n = bn_train_from_file(problem, savefile=big_savefile, do_subset=N, doPrint=True)
			
			bn_save_gmm(gmm_n, gmm_file=gmmfile_n)
			list_gmm.append(gmm_n)
	
	###Now, for each GMM, evaluate Q for each datum in the validation set
	if doPrint:
		print("Evaluating validation samples...",flush=True)
	v = [y+d for y,d in zip(y_samples, d_samples)] #turn into array
	v_traces = [None] * N_val #list of length N_val, each element is a list of length len(N_list)
	for i,vi in enumerate(v):
		if doPrint:
			print("Model evaluation",i+1,'/',len(v),"\t\t",end='\r',flush=True)
		trace = []
		for n_gmm in list_gmm:
			Qvar = gbi_var_of_conditional_pp(n_gmm, vi, verbose=0)#turn into array
			trace.append(Qvar)
		v_traces[i] = trace
	
	###Plot them all, shifting the largest-N estimate of Q to zero to hopefully see many lines converging
	v_traces_adjusted = [[abs(Qvar_n - trace[-1]) for Qvar_n in trace] for trace in v_traces]
	plt.xlabel("N for training GMM")
	plt.ylabel("Difference from final Var[Q] estimate")
	plt.xscale('log')
	for trace in v_traces_adjusted:
		plt.plot(N_list, trace, c='gray', linestyle='dashed')
	plt.show()
	
def bn_compare_model_covariance(problem, datafile, gmmfile, doPrint=True):
	###Load data file
	qoi_train, y_d_train = bn_load_samples(problem, datafile, doPrint)
	#y_d_train = [[yd[0], yd[1]] for yd in y_d_train] #stupid cut for speed
	
	###look at the covariance matrix
	data_together = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
	names_together = ["Q"] + problem.y_names + problem.d_names
	cov = covmatrix_heatmap(data_together, names_together, rescale=True) #gotta rescale to match gmm covariance!
	
	###Load gmm file
	gmm = bn_load_gmm(gmmfile)
	
	###Calculate and plot model covariance
	#https://math.stackexchange.com/questions/195911/calculation-of-the-covariance-of-gaussian-mixtures
	conditional_var_sum = sum([a*C for C,a in zip(gmm.covariances_,gmm.weights_)])
	mu_bar = sum([a * mu for mu,a in zip(gmm.means_,gmm.weights_)])
	var_conditional_mean_sum = sum([a * np.array(mu - mu_bar) * np.array(mu - mu_bar).T for mu,a in zip(gmm.means_,gmm.weights_)])
	total_cov = conditional_var_sum + var_conditional_mean_sum
	
	cov_heatmap(total_cov, names_together)
	
	###Conduct some numerical comparison
	#Frobenius norm: Measures the overall difference between matrices by summing the squared elements of the difference matrix, providing a good general measure of convergence.
	#Operator norm: Measures the maximum linear transformation induced by the difference matrix, useful when interested in the largest potential error in a specific direction.
	#Spectral norm: Measures the largest eigenvalue of the difference matrix, important when considering the impact on principal component analysis.
	#punt it for now
	
def bn_evaluate_model_likelihood(problem, gmmfile, datafile="", samples=[], do_subset=0, doPrint=True):
	if doPrint:
		print("Evaluating model likelihood for",gmmfile,"GMM...",flush=True)
	
	###Get validation set
	if datafile!="" and not samples:
		#use the provided training samples for evaluation
		#this is for determining convergence of the model
		if doPrint:
			print("Loading",datafile,"for training samples...",flush=True)
		qoi_train, y_d_train = bn_load_samples(problem, datafile, doPrint=True, do_subset=do_subset)
		samples = [[q]+yd for q,yd in zip(qoi_train, y_d_train)]
	elif samples and do_subset==0:
		0
	else:
		print("Insufficient instructions for bn_evaluate_model_likelihood")
		sys.exit()
		
	###Load the model to be evaluated
	gmm = bn_load_gmm(gmmfile)
		
	###Standardize the validation set
	gmm_means = gmm.standardized_mean
	gmm_stds = gmm.standardized_std
	standardized_samples = [[(yj-gmm_means[j])/gmm_stds[j] for j,yj in enumerate(y)] for y in samples]
	
	###Calculate the log likelihood of the data
	sample_loglikelihoods = gmm.score_samples(standardized_samples)
	avg_loglikelihood = np.mean(sample_loglikelihoods)

	###print and return
	if doPrint:
		context_str = "validation samples" if datafile=="" else "training set samples"
		print("BN model",gmmfile,"has an average log-likelihood of",avg_loglikelihood,"for",len(standardized_samples),context_str,flush=True)
	return avg_loglikelihood


def bn_measure_likelihood_convergence(problem, big_savefile, doPrint=True):	
	N_list = [1000,4000,7000,10000,40000,70000,100000,400000,700000,1000000]
	
	###Use the savedata to train a series of GMMs, save them if they don't already exist; or, load them if they do
	list_gmm = []
	for N in N_list:
		gmmfile_n = "BN_model_" + str(N) + '.pkl'
		if os.path.exists(gmmfile_n):
			###Add its filename to list
			if doPrint:
				print("Found",gmmfile_n,"...",flush=True)
			list_gmm.append(gmmfile_n)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",gmmfile_n,"...",flush=True)
			gmm_n = bn_train_from_file(problem, savefile=big_savefile, do_subset=N, doPrint=True)
			
			bn_save_gmm(gmm_n, gmm_file=gmmfile_n)
			list_gmm.append(gmmfile_n)
	
	scores = [None] * len(N_list)
	###Now, for each GMM, evaluate the average likelihood of *its* training data
	if doPrint:
		print("Evaluating training data likelihood...",flush=True)
	for i,N in enumerate(N_list):
		#gmmfile: use BN_model_N.pkl
		#datafile: use the big huge one
		#N_val: we're not evaluating validation samples
		#do_subset: N, so that we're evaluating the N gmm with the N values used to train it
		score = bn_evaluate_model_likelihood(problem, gmmfile=list_gmm[i], datafile=big_savefile, do_subset=N, doPrint=True)
		scores[i] = score
	
	###Plot them all, shifting the largest-N estimate of Q to zero to hopefully see many lines converging
	plt.xlabel("N for training GMM")
	plt.ylabel("Average log-likelihood of the model's training data")
	plt.xscale('log')
	plt.plot(N_list, scores, c='r')
	plt.show()
	
def bn_measure_validation_convergence(problem, big_savefile, N_val=0, doPrint=True):	
	N_val = N_val if N_val>0 else 5		
	N_list = [1000,4000,10000,40000,100000,400000,1000000,1635000]
	
	###Get validation set
	#this is for evaluating accuracy of the model
	if doPrint:
		print("Drawing",N_val,"validation samples...",flush=True)
	theta_val = problem.prior_rvs(N_val)
	d_val = problem.sample_d(N_val)
	y_val = [problem.eta(theta, d) for theta,d in zip(theta_val,d_val)]
	q_val = [problem.H(theta) for theta in theta_val]
	samples = [[q_val[i]]+y_val[i]+d_val[i] for i in range(N_val)]
	
	###Use the savedata to train a series of GMMs, save them if they don't already exist; or, load them if they do
	list_gmm = []
	for N in N_list:
		gmmfile_n = "BN_model_" + str(N) + '.pkl'
		if os.path.exists(gmmfile_n):
			###Add its filename to list
			if doPrint:
				print("Found",gmmfile_n,"...",flush=True)
			list_gmm.append(gmmfile_n)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",gmmfile_n,"...",flush=True)
			gmm_n = bn_train_from_file(problem, savefile=big_savefile, do_subset=N, doPrint=True)
			
			bn_save_gmm(gmm_n, gmm_file=gmmfile_n)
			list_gmm.append(gmmfile_n)
	
	scores = [None] * len(N_list)
	###Now, for each GMM, evaluate the average likelihood of *its* training data
	if doPrint:
		print("Evaluating training data likelihood...",flush=True)
	for i,N in enumerate(N_list):
		#gmmfile: use BN_model_N.pkl
		#datafile: we're not evaluating training samples
		#samples: validation set
		#do_subset: we're not evaluating training samples, so its irrelevant
		score = bn_evaluate_model_likelihood(problem, gmmfile=list_gmm[i], samples=samples, doPrint=True)
		scores[i] = score
	
	###Plot them all, shifting the largest-N estimate of Q to zero to hopefully see many lines converging
	plt.xlabel("N for training GMM")
	plt.ylabel("Average log-likelihood of the validation set (Nv="+str()+')')
	plt.xscale('log')
	plt.plot(N_list, scores, c='r')
	plt.show()
	
#Similar to bn_train_from_file
#Except I want to evaluate BIC and MSE, and MAE for each ncomp
def bn_train_evaluate_ncomp(problem, trainfile, valfile, doPrint=True, doPlot=True):
	###Setup
	do_subset=0
	ncomps = [1+ncomp for ncomp in range(49)]
	print("Evaluating",trainfile,"training data for GMM with number of components:",ncomps,flush=True)
	BICs = [None]*len(ncomps)
	scores = [None]*len(ncomps)
	LRTs = [None]*len(ncomps)

	###Load training data file
	qoi_train, y_d_train = bn_load_samples(problem, trainfile, doPrint, do_subset)
	
	###Load validation data file
	qoi_val, y_d_val = bn_load_samples(problem, valfile, doPrint, do_subset)
	val_samples = [[qoi_val[i]]+y_d_val[i] for i in range(len(qoi_val))]
	
	###For each ncomp:
	for i,ncomp in enumerate(ncomps):
		###Train and save a model with the full data with that ncomp
		modelsave = "BN_model_"+str(len(qoi_train))+'_ncomp'+str(ncomp) +'.pkl'
		if os.path.exists(modelsave):
			###Load it
			if doPrint:
				print("Loading",modelsave,"...",flush=True)
			gmm_ncomp = bn_load_gmm(modelsave)
		else:
			###Train and save it
			if doPrint:
				print("Training and saving",modelsave,"...",flush=True)
			gmm_ncomp = gbi_train_model(qoi_train, y_d_train, verbose=2, ncomp=ncomp, careful=True)
			bn_save_gmm(gmm_ncomp, gmm_file=modelsave)
		
		###save the BICs
		p_mean = np.mean(qoi_train)
		p_std = np.std(qoi_train)
		yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_train]
		d_mean = np.mean(y_d_train, axis=0)
		d_std = np.std(y_d_train, axis=0)
		yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_d_train)]
		data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_train)])
		BIC = gmm_ncomp.bic(data)
		BICs[i] = BIC
		print(BIC,flush=True)
		
		###Compare to validation set, grab MSE and MAE
		score = bn_evaluate_model_likelihood(problem, gmmfile=modelsave, samples=val_samples, doPrint=True)
		scores[i] = score

	###go back and calculate the LRTs
	#https://geostatisticslessons.com/lessons/gmm
	for i,score in enumerate(scores):
		if i+1 < len(scores): #this is intended to exclude the final score
			minus2loglambda = 2 * (scores[i+1] - scores[i])
		else: #end of the list
			#gotta get score for g1 = g0 + 1
			modelsave = "BN_model_"+str(len(qoi_train))+'_ncomp'+str(ncomps[-1]+1)+'.pkl'
			if os.path.exists(modelsave):
				###Load it
				if doPrint:
					print("Loading",modelsave,"...",flush=True)
				gmm_g1 = bn_load_gmm(modelsave)
			else:
				###Train and save it
				if doPrint:
					print("Training and saving",modelsave,"...",flush=True)
				gmm_g1 = gbi_train_model(qoi_train, y_d_train, verbose=2, ncomp=ncomps[-1]+1, careful=True)
				bn_save_gmm(gmm_g1, gmm_file=modelsave)
			g1_score = bn_evaluate_model_likelihood(problem, gmmfile=modelsave, samples=val_samples, doPrint=True)
			
			minus2loglambda = 2 * (g1_score - scores[i])
		print(minus2loglambda,flush=True)
		LRTs[i] = minus2loglambda

	if doPlot:
		bn_train_evaluate_ncomp_plot(BICs, LRTs)
		
def bn_train_evaluate_ncomp_plot(BICs, LRTs):
	BICs=[]
	LRTs=[]
	print(BICs,flush=True)
	###Extend with saved data
	BICs.extend([-16691833.485239271, 
	#-156198912.82833186, 
	#-163859555.75851455, 
	#-169702599.09363446, 
	#-184588834.47382435, 
	#-191970333.32038167, 
	#-206948107.51750967, 
	#-205552688.5659083, 
	#-213521927.38937464, 
	-221332299.38899246, 
	#-224094265.49004027, 
	#-224694984.45991743, 
	#-226695631.51102397, 
	#-225555896.76695704, 
	#-228518762.29767773, 
	#-236564864.6959652, 
	#-232971998.46472523,  
	#-233091261.667456, 
	#-235801551.91970456, 
	-236505686.6745928, #20
	#-238863810.6797619, #21
	#-242574931.72879303,
	#-242241982.31413585,
	#-241716697.1986963,
	#-240675759.35915956,
	#-242880225.4465139, #26
	#-244747205.99100965, #27
	#-244505353.69673616, #28
	#-244825887.993535, #29
	-246234260.62625143, #30
	#-246047357.86784902, #31
	#-247755979.07848936, #32
	#-248006778.90827712, #33
	#-246476728.4933763, #34
	#-249868794.03001457, #35
	#-250577132.59812415, #36
	#-249825741.27800956, #37
	#-251868694.3114901, #38
	#-251361870.71319196, #39
	-251633736.07104686, #40
	#-251412903.28971827, #41
	#-254042591.58247107, #42
	#-254065455.92121324, #43
	#-252280005.26191378, #44
	#-251614701.98303768, #45
	#-252855694.75685948, #46
	#-254248535.34322706, #47
	#-254760487.43156114, #48
	#-254193656.49412128, #49
	-256328123.80002046, #50
	
	#-253652561.19815692, #51
	#-254322534.7948813, #52
	#-255520595.29658613, #53
	#-256331530.64211535, #54
	#-256197410.4409594, #55
	
	#-255081159.36240014, #56
	#-256527797.1152506, #57
	#-255586720.07550997, #58
	#-255756605.6004528, #59
	-258218702.1485621, #60
	
	#-257462209.18714023,#61
	#-256818544.98253,#62
	#-257593838.56346256,#63
	#-259556343.32860583,#64
	#-256825583.1222684, #65
	
	#-257497927.4484349,#66
	#-259752913.9879013,#67
	#-260453781.10451147,#68
	#-258016461.481482,#69
	-259624881.49156117, #70
	
	#-257514066.76771572, #71
	#-258676703.9606113, #72
	#-260409247.97088462, #73
	#-258590901.91337222, #74
	#-261407436.56752124, #75
	
	#-261728883.36422002, #76
	#-258718760.19900984, #77
	#-259911816.39266628, #78
	#-261989661.7223139, #79
	-260045696.01722667,#80
	-260782284.5728303 #90
	#100
	-262883685.91923916 #110
	#120
	#130
	])
	"""
	#scores = scores + []
	LRTs.extend([84.64506286831141, #1
	4.527445162701113, #2
	3.7711393767590664, #3
	8.268826614978948, #4
	5.140782007736448, #5
	10.048416146898873, #6
	-1.3586799058341654, #7
	4.764431175551707, #8 
	4.642626800584196, #9
	2.267777451483198, #10
	0.06730020848874574, #11
	1.4827566950605728, #12
	-1.057275376166217, #13
	1.6019359746272244, #14
	4.970770161204086, #15
	-2.2348546755937946, #16
	0.007160972492613382, #17
	1.639094093962683, #18
	0.8380878848860789, #19
	1.2650589385557964, #20
	1.8695235072624428, #21
	0.1912567549178732, #22
	-0.2985982539842098, #23
	-0.7671528181587348, #24
	1.1869351485157722, #25
	1.5097198023039482, #26
	-0.46685502872105644, #27
	0.6694529000964167, #28
	0.663964179556956, #29
	0.2668767063933899, #30
	1.0816540524838274, #31
	0.32215632677261397, #32
	-0.7845671129333596, #33
	1.9800129158032007, #34
	0.3921928244507171, #35
	-0.6801522541597365, #36
	1.544775283556163, #37
	-0.36524014123071424, #38
	0.2538961638286992, #39
	-0.2160321713670612, #40
	1.6032223779436947, #41
	0.20434925449183083, #42
	-1.2555952760467903, #43
	-0.29049415539802226, #44
	0.7927465410666628, #45
	0.7845062146639634, #46
	0.41046871316200395, #47
	-0.40147259034228, #48
	1.3470651331925012, #49
	-1.5768404773600366, #50
	
	0.3795737589875614, #51
	0.794574832903777, #52
	0.5044917697968287, #53
	-0.08964912761081223, #54
	-0.667483869227226, #55
	
	0.9631603213923086, #56
	-0.6959907972261306, #57
	0.05898135158369655, #58
	1.7312252414515399, #59
	-0.7203513954753475, #60
	
	-0.3306525963469653, #61
	0.5986245844209748, #62
	1.4667463496177504, #63
	-1.9345445342086123, #64
	0.4969709854927373, #65
	
	1.481686993224315, #66
	0.34846920754165467, #67
	-1.551469568324677, #68
	0.930682683580585, #69
	-1.1720017312781863, #70
	
	0.8090601572023957, #71
	1.073025183411886, #72
	-1.0998437501094998, #73
	1.6661429504939917, #74
	0.05915066189064078, #75
	
	#76
	#77
	#78
	#79
	#80
	])
	"""
	ncomps = [1,10,20,30,40,50,60,70,80,90,100,110]
	
	###Print, plot
	fig, ax1 = plt.subplots()

	# Plot the first data set on the left y-axis
	ax1.plot(ncomps, BICs, color='blue')
	ax1.set_xlabel("Number of GMM components")
	ax1.set_xticks(ncomps)
	ax1.set_ylabel('BIC', color='blue')
	ax1.tick_params(axis='y', labelcolor='blue')

	# Create the second axes sharing the x-axis
	#ax2 = ax1.twinx()

	# Plot the second data set on the right y-axis
	#ax2.plot(ncomps, LRTs, color='orange')
	#ax2.set_ylabel('LRT', color='orange')
	#ax2.tick_params(axis='y', labelcolor='orange')
	plt.show()

def bn_train_evaluate_ncomp_sanitycheck(problem, trainfile, valfile, doPrint=True, doPlot=True):
	###Setup
	do_subset=0
	ncomps = [1,10,20,30,40,50,60,70,80,90]
	print("Evaluating",trainfile,"training data for GMM with number of components:",ncomps,flush=True)
	BICs = [None]*len(ncomps)
	scores = [None]*len(ncomps)
	penalties = [None]*len(ncomps)

	###Load training data file
	qoi_train, y_d_train = bn_load_samples(problem, trainfile, doPrint, do_subset)
	
	###Load validation data file
	qoi_val, y_d_val = bn_load_samples(problem, valfile, doPrint, do_subset)
	val_samples = [[qoi_val[i]]+y_d_val[i] for i in range(len(qoi_val))]
	
	###For each ncomp:
	for i,ncomp in enumerate(ncomps):
		###Train and save a model with the full data with that ncomp
		modelsave = "BN_model_"+str(len(qoi_train))+'_ncomp'+str(ncomp) +'.pkl'
		if os.path.exists(modelsave):
			###Load it
			if doPrint:
				print("Loading",modelsave,"...",flush=True)
			gmm_ncomp = bn_load_gmm(modelsave)
		else:
			###Train and save it
			print("Model",modelsave,"doesn't exist.",flush=True)
			sys.exit()
		
		###save the BICs - VAL
		p_mean = np.mean(qoi_val)
		p_std = np.std(qoi_val)
		yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_val]
		d_mean = np.mean(y_d_val, axis=0)
		d_std = np.std(y_d_val, axis=0)
		yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_d_val)]
		data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_val)])
		
		BIC = gmm_ncomp.bic(data)
		BICs[i] = BIC
		score = -2 * gmm_ncomp.score(data) * data.shape[0]
		scores[i] = score
		penalty = gmm_ncomp._n_parameters() * np.log(data.shape[0])
		penalties[i] = penalty
		print(BIC,gmm_ncomp.score(data),gmm_ncomp._n_parameters(),flush=True)
		#_n_parameters calculates K*((D*D - D)/2 + 2D + 1)-1, for D=62

	if doPlot:
		fig, ax1 = plt.subplots()

		# Plot the first data set on the left y-axis
		ax1.plot(ncomps, BICs, color='blue')
		ax1.set_xlabel("Number of GMM components")
		ax1.set_xticks(ncomps)
		ax1.set_ylabel('BIC', color='blue')
		ax1.tick_params(axis='y', labelcolor='blue')
		ax1.plot(ncomps, penalties, color='red')
		ax1.plot(ncomps, scores, color='green')
		plt.show()




#Do bootstrap sampling & evaluation to find 95% confidence intervals
def bn_train_convergence_confidence(problem, trainfile, N_bootstrap, ncomp, startnum=0, doPrint=True):
	###Load training data file
	qoi_train, y_d_train = bn_load_samples(problem, trainfile, doPrint, do_subset)
	
	###Iterate
	i_list = range(startnum, N_bootstrap)
	for i in i_list:
		###Make bootstrap sample
		i_samples = np.random.randint(0, len(qoi_train), size=len(qoi_train))
		qoi_bootstrap = [qoi_train[j] for j in i_samples]
		y_d_bootstrap = [y_d_val[j] for j in i_samples]
		
		###Train and save GMM on that
		modelsave = "BN_model_"+str(len(qoi_train))+'_bootstrap'+str(i)+'.pkl'
		if os.path.exists(modelsave):
			print("GMM",modelsave,"already exists, terminating.")
			sys.exit()
		else:
			if doPrint:
				print("Training and saving",modelsave,"...",flush=True)
			gmm_i = gbi_train_model(qoi_bootstrap, y_d_bootstrap, verbose=2, ncomp=ncomp, careful=True)
			bn_save_gmm(gmm_i, gmm_file=modelsave)
		
def bn_evaluate_convergence_confidence(problem, valfile, N_bootstrap, doPrint=True):
	###Load validation data file
	qoi_val, y_d_val = bn_load_samples(problem, valfile, doPrint, do_subset)

	p_mean = np.mean(qoi_bootstrap)
	p_std = np.std(qoi_bootstrap)
	yp_sample = [(qoi - p_mean)/p_std for qoi in qoi_bootstrap]
	d_mean = np.mean(y_d_bootstrap, axis=0)
	d_std = np.std(y_d_bootstrap, axis=0)
	yd_sample = [list((yd - d_mean)/d_std) for i,yd in enumerate(y_d_bootstrap)]
	val_data = np.array([np.hstack([yp_sample[i],yd_sample[i]]) for i,_ in enumerate(qoi_bootstrap)])

	###Iterate
	bootstrapped_scores = [None]*N_bootstrap
	for i in range(N_bootstrap):
		###Load GMM
		modelsave = "BN_model_"+str(len(qoi_train))+'_bootstrap'+str(i)+'.pkl'
		if os.path.exists(modelsave):
			###Load it
			if doPrint:
				print("Loading",modelsave,"...",flush=True)
			gmm_g1 = bn_load_gmm(modelsave)
		else:
			print("error")
			sys.exit()
		
		###Calculate score (average likelihood of validation set), save
		gmm_i.score(val_data)
	
	##Do statistics on the scores
	#Mean
	score_mean = np.mean(bootstrapped_scores)
	
	#stddev
	score_stddev = np.std(bootstrapped_scores)
	
	#conf interval
	def conf_interval(data, conf_level):
		#returns a tuple of the low and high bounds
		return scipy.stats.t.interval(confidence=conf_level, df=len(data)-1, loc=np.mean(data), scale=scipy.stats.sem(data))
		
	conf_intervals = conf_interval(bootstrapped_scores,0.95)
	conf_interval_width = conf_intervals[1] - conf_intervals[0]
	
	return score_mean, score_stddev, conf_interval_width

