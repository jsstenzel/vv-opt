import sys
import os
import scipy.stats
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
import csv
import dill
matplotlib.use('TkAgg')

sys.path.append('.')
#focal plane
from problems.fp_verification_split.fp_problem import *
#analysis
from obed.mcmc import *
from obed.obed_huan import *
from uq.uncertainty_propagation import *
#optimization
from opt.bayesian_opt import *
from approx.chaospy_pce import *

################################
#Useful definitions
################################

theta_nominal = [1.1, 2.5, .001]
QoI_nominal = fp.H(theta_nominal)

d_historical = [
				20,   #t_gain
				30,   #I_gain
				1,	#n_meas_rn
				8,	#d_num
				9600, #d_max
				2	 #d_pow   #approx
			   ]
			   
d_best = [
				600,   #t_gain
				100,   #I_gain
				50,	#n_meas_rn
				100,	#d_num
				12000, #d_max
				3	 #d_pow   #approx
		]
			   
d_worst = [
				1,   #t_gain
				1,   #I_gain
				1,	#n_meas_rn
				2,	#d_num
				1, #d_max
				0.1	 #d_pow   #approx
		]
		
y_nominal = fp._internal_G(dict(zip(fp.theta_names, theta_nominal)), dict(zip(fp.d_names, d_historical)), dict(zip(fp.x_names, fp.x_default)))
#print(y_nominal)

################################
#Analysis functions
################################

def recursive_list_print(input, level=0):
	indent = "\t"*level
	#if type(input) is list:
	if hasattr(input, '__iter__'):
		print(indent+"list["+str(len(input))+"]:")
		for item in input:
			recursive_list_print(item, level=level+1)
	else:
		print(indent+str(input))

def study0_compare_nolikelihood_methods():
	#do direct model evaluations of eta, no PCE
	#Fix total number of eta evaluations at 60,000
	
	##SUBSET 1: model re-use in the 2nd loop
	#direct, likelihood kde, kde over y
	#N = = n_out + n_out*n_3rd
	n_out = 240
	n_3rd = 10000#249
	u,_ = U_Huan_direct_3loop(d_historical, fp, n_out=n_out, n_3rd=n_3rd, n_in=0, doPrint=True)
	print(u)

	#direct, likelihood kde, kde over ytheta, no kde sample reuse
	#N = n_out + n_3rd
	n_out = 30000
	n_3rd = 30000
	u,_ = U_Huan_direct_parallelloop(d_historical, fp, n_out=n_out, n_3rd=n_3rd, n_in=0, doPrint=True)
	print(u)
	#best?
	
	#direct, likelihoo kde, kde over ytheta, kde sample reuse
	#N = n_out
	n_out = 60000
	u,_ = U_Huan_direct_parallelloop(d_historical, fp, n_out=n_out, n_in=0, n_3rd=0, doPrint=True)
	print(u)
	#best?

	##SUBSET 1: no model re-use in the 2nd loop
	#direct, likelihood kde, kde over y
	#N = n_out + n_out*n_in + n_out*n_3rd
	n_out = 50
	n_3rd = 10000#40
	n_in = 39
	u,_ = U_Huan_direct_3loop(d_historical, fp, n_out=n_out, n_3rd=n_3rd, n_in=n_in, doPrint=True)
	print(u)
	#this will obviously be bad

	#direct, likelihood kde, kde over ytheta, no kde sample reuse
	#N = n_out + n_3rd + n_out*n_in
	n_out = 200
	n_3rd = 20000
	n_in = 100
	u,_ = U_Huan_direct_parallelloop(d_historical, fp, n_out=n_out, n_in=n_in, n_3rd=n_3rd, doPrint=True)
	print(u)
	
	#direct, likelihood kde, kde over ytheta, kde sample reuse
	#N = n_out + n_out*n_in
	n_out = 240
	n_in = 249
	u,_ = U_Huan_direct_parallelloop(d_historical, fp, n_out=n_out, n_in=n_in, n_3rd=0, doPrint=True)
	print(u)

	###See if I can rule out a few of those as being way to high error, then consider a subset for convergence study
	###and also varying the nin, nout,n3rd?
	
def study1_u_check(n_out):
	u,ulist,likelihoods,evidences,skipcount = U_Huan(d_historical, fp, fp.G, n_out=n_out, n_in=100, doPrint=True)
	print(u)
	print(likelihoods)
	print(evidences)
	print(skipcount)

"""
look at convergence by comparing like to like n_in
should be able to see the bias?
"""
def study1_convergence(nlist=[]):
	#nlist = [100,300,500,1000,3000]
	ulist_reuse = [4.4749838453835835,5.5030513985692355,5.844435063165167,6.385123526117919,6.909775414327681,7.089762391395981]
	ulist_resample = [15.52263587153617, 11.785538698774506, 11.49525551448308, 9.700575983993778, 8.45884769782254,8.024641932921616]
	skiplist=[]
	for n in nlist:
		u,_,_,_,_ = U_Huan(d_historical, fp, fp.G, n_out=n, n_in=0, doPrint=True)
		print(u)
		ulist_reuse.append(u)
	
	for n in nlist:
		u,_,_,_,skipcount = U_Huan(d_historical, fp, fp.G, n_out=n, n_in=n, doPrint=True)
		print(u)
		ulist_resample.append(u)
		skiplist.append(skipcount)
		
	print(skiplist)
	plt.plot(nlist,ulist_reuse,c='b')
	#[n-skip for n,skip in zip(nlist,skiplist)]
	plt.plot(nlist, ulist_resample,c='r')
	plt.xlabel("n_out")
	plt.ylabel("estimated U")
	plt.show()
	
def study1_convergence_rescale_plot():
	ulist_reuse = [4.4749838453835835,5.5030513985692355,5.844435063165167,6.385123526117919,6.909775414327681,7.089762391395981]
	nlist = [100,300,500,1000,3000,6000]
	
	ulist_resample = [15.52263587153617, 11.785538698774506, 11.49525551448308, 9.700575983993778, 8.45884769782254,8.024641932921616]
	n2list = [n**2 for n in nlist]
	
	plt.plot(nlist, ulist_reuse,c='b')
	plt.plot(n2list, ulist_resample,c='r')
	plt.xlabel("total model evals")
	plt.ylabel("estimated U")
	plt.xscale('log')
	plt.show()
	
	
"""
replicate Appendix B here, varying n_in for fixed n_out
"""
def study1_convergence_deviation(ntries, n_out, nlist = [100,300,500,1000,3000]):
	ulist_reuse = []
	ulist_resample = []

	#for n in nlist:
	#	ulist_i = []
	#	for j in range(ntries):
	#		u,_,_,_,_ = U_Huan(d_historical, fp, fp.G, n_out=n_out, n_in=n, doPrint=False, doReuse=True)
	#		ulist_i.append(u)
	#		print(str(j+1)+"/"+str(ntries),str(u)+'\t', flush=True, end='\r')
	#	print(ulist_i,flush=True)
	#	ulist_reuse.append(ulist_i)
	
	for n in nlist:
		ulist_i = []
		for j in range(ntries):
			u,_,_,_,skip = U_Huan(d_historical, fp, fp.G, n_out=n_out, n_in=n, doPrint=False)
			ulist_i.append(u)
		print(ulist_i,flush=True)
		ulist_resample.append(ulist_i)
		
	ulist_reuse_avg = [np.mean(ui) for ui in ulist_reuse]
	ulist_resample_avg = [np.mean(ui) for ui in ulist_resample]
	ulist_reuse_std = [np.std(ui, ddof=1) for ui in ulist_reuse]
	ulist_resample_std = [np.std(ui, ddof=1) for ui in ulist_resample]
	print(ulist_reuse_std)
	print(ulist_resample_std)
		
	#plt.plot(nlist,ulist_reuse_avg,c='b')
	plt.errorbar(nlist,ulist_reuse_avg,yerr=ulist_reuse_std, ls='--', capsize=4, color='b', ecolor="b")
	#plt.plot(nlist,ulist_resample_avg,c='r')
	plt.errorbar(nlist,ulist_resample_avg,yerr=ulist_resample_std, ls='--', capsize=4, color='r', ecolor="r")
	#plt.savefig("longrun.png")
	plt.show()
	
def study1_design_space_check():
	ubest,_,_,_,_ = U_Huan(d_best, fp, fp.G, n_out=100, n_in=5000, doReuse=True, doPrint=True)
	uworst,_,_,_,_ = U_Huan(d_worst, fp, fp.G, n_out=100, n_in=5000, doReuse=True, doPrint=True)
	print(ubest)
	print(uworst)
	
def study1_optimization(n_init_opt=50, n_iter_opt=100, n_out=2000):
	#define parameters to ensure both algorithms are equal?
	#this will test whether we can accept the bias of reuse for a given
	n_out_reuse = n_out
	n_out_resample = n_out
	n_in_resample = n_out

	#run optimizer for both algorithms
	def U_reuse(d):
		u,_,_,_,_ = U_Huan(d, fp, fp.G, n_out=500, n_in=3000, doReuse=True, doPrint=True)
		if math.isnan(u) or math.isinf(u):
			return 0
		else:
			return u
		
	def U_resample(d):
		u,_,_,_,_ = U_Huan(d, fp, fp.G, n_out=n_out_resample, n_in=n_in_resample, doReuse=False, doPrint=False)
		if math.isnan(u) or math.isinf(u):
			return 0
		else:
			return u
	
	dopt_reuse, Uopt_reuse = bayesian_opt(
		fp, U_reuse, n_init=n_init_opt, n_iter=n_iter_opt, f_noise="gaussian", findMinimum=False, plotConv=True)
	
	#dopt_resample, Uopt_resample = bayesian_opt(
	#	fp, U_resample, n_init=n_init_opt, n_iter=n_iter_opt, f_noise="gaussian", findMinimum=False, plotConv=True)
	
	#hypothesis: we get the same answer for both - proving that reuse is more efficient for the same result
	#note that for maximization, we shouldnt care about
	print("Optimal d under reuse:",dopt_reuse)
	print("\tHighest utility:", Uopt_reuse)
	#print("Optimal d under resample:",dopt_resample)
	#print("\tHighest utility:", Uopt_resample)

def study15_convergence_randomizer(n_out=10000):
	nlist = [100,300,500,1000]#,3000]
	ulist_resample = [14.505063968535362,12.233080207486598,10.83320792,9.374276209]
	ulist_rerandom = [14.91182393,11.94412013,10.74615738,9.356828537]
	#for n in nlist:
	#	u,_,_,_,skipcount = U_Huan(d_historical, fp, fp.G, n_out=n_out, n_in=n, doReuse=True, doPrint=True)
	#	print(u)
	#	ulist_resample.append(u)
		
	#for n in nlist:
	#	u,_,_,_,skipcount = U_Huan(d_historical, fp, fp.G, n_out=n_out, n_in=n, doReuse=2, doPrint=True)
	#	print(u)
	#	ulist_rerandom.append(u)
		
	plt.plot(nlist,ulist_resample,c='r')
	plt.plot(nlist, ulist_rerandom,c='orange')
	plt.xlabel("n_in")
	plt.ylabel("estimated U")
	plt.show()


def study2_likelihood_pce_check():
	np.set_printoptions(suppress=True)
	#First, we need to develop the PCE model estimate, and the associated likelihood fns
	
	#now calculate the utility for an example design, look at convergence
	approx_lin = likelihood_G_PCE_linregress(fp, n_eval=1000, order=3)
	approx_psp, jointdist = likelihood_G_PCE_psp(fp, 6, 8)
	print(len(approx_lin))
	print(len(approx_psp))

	thetad = theta_nominal + d_historical
	approx_evals = [poly(*thetad) for poly in approx_psp]
	print(approx_evals)
	print("PSP sum:",np.sum(approx_evals,axis=0)) #psp
	print("PSP avg:",np.mean(approx_evals,axis=0)) #psp
	print("Linear regression model:",approx_lin(*thetad)) #lin
	print("True model",fp.G(theta_nominal,d_historical)) #true

"""
to measure error here, I should develop the linear PCE model for a variety of different ni's
for each ni, I should evaluate by taking 10,000 samples of the d and theta space, evaluate with both models, look at the distribution of the difference,
and finally, plot the average difference for each ni
"""
def study2_likelihood_pce_convergence(nlist, test_N=1000):
	avg_diffs = []
	for ni in nlist:
		print(ni,flush=True)
		#construct
		approx_lin = likelihood_G_PCE_linregress(fp, n_eval=ni, order=3)
		
		#run the same test_N samples through both model and original
		thetas = fp.prior_rvs(test_N)
		ds = fp.sample_d(test_N)
		thetad_list = [np.concatenate([t,d]) for t,d in zip(thetas,ds)]
		approx_evals = [approx_lin(*thetad) for thetad in thetad_list]
		direct_evals = [fp.G(t,d) for t,d in zip(thetas,ds)]
		
		diffs = [[abs(approx_i[j] - direct_i[j]) for approx_i,direct_i in zip(approx_evals,direct_evals)] for j in range(len(direct_evals[0]))]
		avg_diff = [np.mean(diff_i) for diff_i in diffs]
		print(avg_diff,flush=True)
		avg_diffs.append(avg_diff)
	
	print(avg_diffs)

	plt.plot(nlist,[avg[0] for avg in avg_diffs],c='g')
	plt.yscale('log')
	plt.show()
	plt.plot(nlist,[avg[1] for avg in avg_diffs],c='g')
	plt.yscale('log')
	plt.show()
	plt.plot(nlist,[avg[2] for avg in avg_diffs],c='g')
	plt.yscale('log')
	plt.show()

	
def study2_convergence_deviation(ntries=30):
	nlist = [1000]
	ulist_direct = []
	ulist_pce = []
	
	#for ni in nlist:
	#	ulist_i = []
	#	for j in range(ntries):
	#		u,_,_,_,_ = U_Huan(d_historical, fp, fp.G, n_out=ni, n_in=0, doPrint=False, doReuse=True)
	#		ulist_i.append(u)
	#		print(str(j+1)+"/"+str(ntries),str(u)+'\t', flush=True, end='\r')
	#	print(ulist_i,flush=True)
	#	ulist_direct.append(ulist_i)
		
	approx_lin = likelihood_G_PCE_linregress(fp, n_eval=2000, order=3) #this is based on the study2_likelihood_pce_convergence test
	def approx_G(theta,d):
		thetad = np.concatenate([theta,d])
		return approx_lin(*thetad)
	
	for ni in nlist:
		ulist_i = []
		for j in range(ntries):
			u,_,_,_,skip = U_Huan(d_historical, fp, approx_G, n_out=ni, n_in=0, doPrint=True, doReuse=True)
			ulist_i.append(u)
		print(ulist_i,flush=True)
		ulist_pce.append(ulist_i)
		
	ulist_direct_avg = [np.mean(ui) for ui in ulist_direct]
	ulist_pce_avg = [np.mean(ui) for ui in ulist_pce]
	ulist_direct_std = [np.std(ui, ddof=1) for ui in ulist_direct]
	ulist_pce_std = [np.std(ui, ddof=1) for ui in ulist_pce]
	print(ulist_direct_std)
	print(ulist_pce_std)
		
	#plt.plot(nlist,ulist_reuse_avg,c='b')
	plt.errorbar(nlist,ulist_direct_avg,yerr=ulist_direct_std, ls='--', capsize=4, color='r', ecolor="r")
	#plt.plot(nlist,ulist_resample_avg,c='r')
	plt.errorbar(nlist,ulist_pce_avg,yerr=ulist_pce_std, ls='--', capsize=4, color='g', ecolor="g")
	plt.show()


"""
make a nlist
assume we're doing model reuse, and we have some N=1000 for our obed monte carlo
i want to investigate the marginal utility of spending more model evals on PCE, as opposed to monte carlo
for each ni, calculate:
the uncertainty of U for N+ni steps
the uncertainty of U for N steps, with a PCE trained on ni model evals
product will be two plots showing U converging at different rates??
"""
def study2_compare_pce_Uvar(N, nlist):
	u_directs = []
	u_pces = []
	
	u_N,_ = U_Huan(d_historical, fp, fp.G, n_out=N, n_in=0, doPrint=True)
	u_directs.append(u_N)
	
	for ni in nlist:
		#direct
		u_direct,_ = U_Huan(d_historical, fp, fp.G, n_out=N+ni, n_in=0, doPrint=True)
		u_directs.append(u_direct)
	
		#pce
		approx_lin = likelihood_G_PCE_linregress(fp, n_eval=ni, order=3)
		def approx_G(theta,d):
			thetad = np.concatenate([theta,d])
			return approx_lin(*thetad)
			
		u_pce,_ = U_Huan(d_historical, fp, approx_G, n_out=N, n_in=0, doPrint=True)
		u_pces.append(u_pce)
		
	#plot
	ns = [N+ni for ni in nlist]
	plt.plot([N]+ns,u_directs,c='r')
	plt.plot(ns,u_pces,c='g')
	plt.show()
	
import time
	
def study2_compare_time():
	np.set_printoptions(suppress=True)
	#First, we need to develop the PCE model estimate, and the associated likelihood fns
	thetad = theta_nominal + d_historical
		
	t_0 = time.perf_counter()
	
	#now calculate the utility for an example design, look at convergence
	approx_lin = likelihood_G_PCE_linregress(fp, n_eval=1000, order=1)
	#print(approx_lin)
	
	t_1 = time.perf_counter()

	#PCE model
	for _ in range(1000):
		approx_lin(*thetad) #lin
	
	t_2 = time.perf_counter()
	
	#true model
	for _ in range(1000):
		fp.G(theta_nominal,d_historical) #true
	
	t_3 = time.perf_counter()
	
	print("initial PCE time:",t_1 - t_0)
	print("time for 1000 PCE evals:",t_2 - t_1)
	print("time for 1000 direct evals:",t_3 - t_2)
	
	
def study2_optimization(n_init_opt=50, n_iter_opt=100, n_out=2000):
	#define parameters to ensure both algorithms are equal?
	#this will test whether we can accept the bias of reuse for a given
	n_out_reuse = n_out
	n_out_resample = n_out
	n_in_resample = n_out

	#run optimizer for both algorithms
	def U_direct(d):
		u,_,_,_,_ = U_Huan(d, fp, fp.G, n_out=n_out_reuse, n_in=0, doReuse=True, doPrint=False)
		return u
	
	approx_lin = likelihood_G_PCE_linregress(fp, n_eval=ni, order=3)
	def approx_G(theta,d):
		thetad = np.concatenate([theta,d])
		return approx_lin(*thetad)	

	def U_pce(d):
		u_pce,_,_,_,_ = U_Huan(d, fp, approx_G, n_out=n_out, n_in=0, doReuse=True, doPrint=True)
		return u_pce
	
	dopt_direct, Uopt_direct = bayesian_opt(
		fp, U_direct, n_init=n_init_opt, n_iter=n_iter_opt, f_noise="gaussian", findMinimum=False, plotConv=True)
	
	dopt_pce, Uopt_pce = bayesian_opt(
		fp, U_pce, n_init=n_init_opt, n_iter=n_iter_opt, f_noise="gaussian", findMinimum=False, plotConv=True)
	
	#hypothesis: we get the same answer for both - proving that reuse is more efficient for the same result
	#note that for maximization, we shouldnt care about
	print("Optimal d for direct model eval:",dopt_direct)
	print("\tHighest utility:", Uopt_direct)
	print("Optimal d for PCE model approx:",dopt_pce)
	print("\tHighest utility:", Uopt_pce)


if __name__ == '__main__':  
	#SETUP
	
	
	###STUDY 0: think about what i can do if i can't evaluate the likelihood
	#study0_compare_nolikelihood_methods()
	"""
	###STUDY 1
	study1_u_check(1000)
	
	study1_convergence([6000])
	study1_convergence_rescale_plot()
	
	study1_convergence_deviation(ntries=50, n_out=1000)	
	
	study1_design_space_check()
	study1_optimization()
	#run optimizer for both algorithms
	#hypothesis: we get the same answer for both - proving that reuse is more efficient for the same result
	
	
	###STUDY 1.5 - does random slicing of n_out for reuse decrease bias?
	study15_convergence_randomizer()
	
	
	###STUDY 2 -- then compare those to pce?
	study2_likelihood_pce_check()
	#ill just go with the linar PCE
	
	#the thing i want to show is the gains in U_predict you get by putting n model evaluations to the PCE solver, as opposed to the simulation
	study2_likelihood_pce_convergence(nlist=[10,50,100,500,1000,5000])
	
	#study2_convergence_deviation() #didnt plot this

	study2_compare_pce_Uvar(100, nlist=[100,200,500,750,1000,2000])
	
	study2_compare_time()
	
	#study2_optimization() #didnt run this
	#run optimizer, with reuse, for both direct and pce G
	#hypothesis: we get the same answer for both - proving that pce is more efficient for the same result
	"""
	study1_convergence_deviation(ntries=50, n_out=5000, nlist = [5000])



"""
240/240 [-1.91763571]
-1.9797318251454674
30000/30000 [3.13380215]]]]
2.70271473205692
60000/60000 [3.12252885]]5]
3.0006725455507386
50/50 [-3.27811644]
-2.015610985865101
200/200 [0.60322796]]
2.5660640333031015
240/240 [3.54941106]
3.3703073224466853
"""