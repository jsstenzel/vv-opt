import numpy as np
import scipy.stats

#https://quant.stackexchange.com/questions/21764/stopping-monte-carlo-simulation-once-certain-convergence-level-is-reached

#determine the number of iterations of algorithm necessary to find convergence in the test_statistic
#keep going until 
#return n
def solve_convergence(algorithm, ci=0.95, closeness=0.95, min_runs=100, **kwargs):
	q = scipy.stats.norm.ppf(ci) #phi^-1(ci), quantile fn of x. for example, ci=0.95 means q=1.96
	threshold = 1.0-closeness #this is how many mean-units we want our CI to be
	
	#Start by running for the minimum number of runs
	N = min_runs
	init_sample = []
	for _ in range(N):
		init_sample.append(algorithm(kwargs))
	
	#Start with the initial sample mean and variance
	iter_sample_mean = np.mean(init_sample)
	iter_sample_var = np.var(init_sample, ddof=1) #its a sample variance
	
	while q*np.sqrt(iter_sample_var/N) > threshold*iter_sample_mean:
		new_val = algorithm(kwargs)
		
		prev_mean = iter_sample_mean
		iter_sample_mean = prev_mean + (new_val - prev_mean)/(N+1)
		iter_sample_var = (1-(1/N))*iter_sample_var + (n+1)*(iter_sample_mean - prev_mean)**2
		
		N += 1
		
	return N, iter_sample_mean
	

#This is just like the above, but can be used in a convergence solver for a more complicated case
def is_algorithm_converged(sample, ci=0.95, closeness=0.95):
	q = scipy.stats.norm.ppf(ci) #phi^-1(ci), quantile fn of x. for example, ci=0.95 means q=1.96
	threshold = 1.0-closeness #this is how many mean-units we want our CI to be
	N = size(sample)
	
	#Find sample mean and variance of the data
	sample_mean = np.mean(sample)
	sample_var = np.var(sample, ddof=1) #its a sample variance
	
	#without knowing anything else, report whether this data meets the convergence criterion
	return q*np.sqrt(sample_var/N) <= threshold*sample_mean
	
	
#determine the number of iterations of algorithm necessary to find convergence in the test_statistic
#keep going until a bootstrap finds that ci around the statistic estimate x is x +- x*(1-closeness)
#return n
def solve_bootstrapped_convergence(algorithm, ci=0.95, closeness=0.95, **kwargs):
	#at some point, make a similar algorithm to the above and try doing bootstrap estimate of var instead?