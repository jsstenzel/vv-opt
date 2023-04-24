import numpy as np
import scipy.stats


#####################################################################
# Generic algorithms for observing convergence of the 3 subproblems
#####################################################################

#To try to answer how large n3 needs to be, im going to try plotting trace of gaussian kde eval'd & summed at the same points
#This assumes we're using scipy kde solver like in calc_likelihood_kernel -- i probably want to look at different algs
def trace_scipy_kde(d, exp_fn, H_fn, G_fn, p_theta_rvs, p_theta_pdf, nmax=10000):
    #nmax = 10000
    #d = 0.5
    thetas = p_theta_rvs(100)
    ys = [exp_fn(theta, d) for theta in thetas]
    
    test_thetas = p_theta_rvs(100)
    test_y_theta = [[exp_fn(theta, d),theta] for theta in test_thetas]
    trace = []
    for n in range(nmax):
        thetas = np.append(thetas, p_theta_rvs(1)[0])
        ys = np.append(ys, eta(thetas[-1], d))
        
        y_theta_values = np.vstack([ys, thetas])
        likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values)
        
        #grab metrics from this kde
        #i want to trace some kind of change from one step to the next
        test_calcs = [likelihood_kernel.pdf(pair)[0] for pair in test_y_theta]
        
        trace_metric = sum(test_calcs)
        trace.append(trace_metric)
    #print(trace, flush=True)
    plt.plot(trace)
    plt.show()
	
	
#How many MCMC steps are needed to generate a good sample of the posterior?
#This function doesnt seem to work :/
#assumes u = var(H(posterior))
def loop2_convergence(d, y, exp_fn, H_fn, G_fn, p_theta_rvs, p_theta_pdf, burnin=0, lag=1, min=100, doLog=True, doPlot=True, doPrint=False):
	#basically i want to go until Ui is converged

	#p_theta_rvs prior function
	#p_theta_pdf prior function
	#exp_fn experiment model / eta
	#H_fn high-level system mdel / QoI function
	#d = 0.5
	#y = 0.6 just pick one from likelihood, or use a plausible fixed value to have comparable results
	#min the number of mcmc steps to start with when calculating the utility v/ H variance
	
	if not doLog and not doPlot and not doPrint:
		print("No, you probably want loop2_convergence to either log or plot or print. Otherwise that would be a waste!", flush=True)
		exit()
	else:
		notify = "Doing" + ("", " logging,")[doLog] + ("", " plotting,")[doPlot] + ("", " print all,")[doPrint]
		print(notify[:-1])
	
	#assume we have a good likelihood kernel to use in loop3. set it up to reuse here
	n1 = 10000
	likelihood_kernel, _ = calc_likelihood_kernel(d, exp_fn, p_theta_rvs, n1, showplot=False)
	
	#setup loop
	theta_current = p_theta_rvs(1)[0]
	mcmc_trace = []
	U_trace = []
	i = 0
	
	while not is_algorithm_converged(U_trace, min, ci=0.975, closeness=0.99):
		#one MCMC loop step
		#Propose a new value of theta from the prior
		theta_proposal = p_theta_rvs(1)[0]
		#Compute likelihood of the "data" for both thetas - reusing the loop-1 likelihood evaluations
		likelihood_current = likelihood_kernel.evaluate([y,theta_current]) #using kde to estimate pdf
		likelihood_proposal = likelihood_kernel.evaluate([y,theta_proposal]) #calc probability at the data y
		#Compute acceptance ratio
		prior_current = p_theta_pdf(theta_current)
		prior_proposal = p_theta_pdf(theta_proposal)
		p_current = likelihood_current * prior_current
		p_proposal = likelihood_proposal * prior_proposal
		R = p_proposal / p_current
		#Accept our new value of theta according to the acceptance probability R:
		if np.random.random_sample() < R:
			theta_current = theta_proposal
		#Include theta_current in my trace according to rules
		i += 1
		if not( i > burnin and i%lag == 0):
			continue #i.e. run the mcmc part again, dont add corresponding U to the trace
		mcmc_trace.append(theta_current)
		
		#We also want to use min here to make sure that we're not taking a meaningless variance for our utility, with too little mcmc data
		if len(mcmc_trace) <= min:
			continue
			
		#Now, use my posterior samples to calculate H(theta|y,d) samples
		H_theta_posterior = [H_fn(theta) for theta in mcmc_trace]
		
		#Now find the variance of that. Return the inverse as the utility; consider cost separately
		var_H = np.var(H_theta_posterior, ddof=1) #its a sample variance
		U = (1.0/var_H)
		U_trace.append(U)
		#print(len(U_trace), end='\r', flush=True)
		
		#more informative print statement
		if doPrint:
			sample_mean = np.mean(U_trace)
			sample_var = np.var(U_trace, ddof=1)
			N_latest = len(U_trace)
			print(N_latest, U, sample_mean, sample_var, scipy.stats.norm.ppf(0.95), scipy.stats.norm.ppf(0.95)*np.sqrt(sample_var/N_latest)/sample_mean, flush=True)
		
	print("N for convergence is", len(U_trace), ", mean U is", np.mean(U_trace))
	
	
#To try to answer how large n2 needs to be, im going to try plotting trace of mcmc and resulting utility
#assumes u = var(H(posterior))
def trace_loop2(d, y, exp_fn, H_fn, G_fn, p_theta_rvs, p_theta_pdf, burnin=0, lag=1, min=100, doLog=True, doPlot=True):
	#p_theta_rvs prior function
	#p_theta_pdf prior function
	#exp_fn experiment model / eta
	#H_fn high-level system mdel / QoI function
	#d = 0.5
	#y = 0.6 just pick one from likelihood, or use a plausible fixed value to have comparable results
	#min the number of mcmc steps to start with when calculating the utility v/ H variance
	
	if not doLog and not doPlot:
		print("No, you probably want trace_loop2 to either log or plot. Otherwise that would be a waste!", flush=True)
		exit()
	
	#assume we have a good likelihood kernel to use in loop3. set it up to reuse here
	n1 = 10000
	likelihood_kernel, _ = calc_likelihood_kernel(d, exp_fn, p_theta_rvs, n1, showplot=False)
	
	#setup loop
	theta_current = p_theta_rvs(1)[0]
	mcmc_trace = []
	U_trace = []
	i = 0
	
	n2 = 100000
	for _ in range(n2):
		#one MCMC loop step
		#Propose a new value of theta from the prior
		theta_proposal = p_theta_rvs(1)[0]
		#Compute likelihood of the "data" for both thetas - reusing the loop-1 likelihood evaluations
		likelihood_current = likelihood_kernel.evaluate([y,theta_current]) #using kde to estimate pdf
		likelihood_proposal = likelihood_kernel.evaluate([y,theta_proposal]) #calc probability at the data y
		#Compute acceptance ratio
		prior_current = p_theta_pdf(theta_current)
		prior_proposal = p_theta_pdf(theta_proposal)
		p_current = likelihood_current * prior_current
		p_proposal = likelihood_proposal * prior_proposal
		R = p_proposal / p_current
		#Accept our new value of theta according to the acceptance probability R:
		if np.random.random_sample() < R:
			theta_current = theta_proposal
		#Include theta_current in my trace according to rules
		i += 1
		if not( i > burnin and i%lag == 0):
			continue #i.e. run the mcmc part again, dont add corresponding U to the trace
		mcmc_trace.append(theta_current)
		
		#We also want to use min here to make sure that we're not taking a meaningless variance for our utility, with too little mcmc data
		if len(mcmc_trace) <= min:
			continue
			
		#Now, use my posterior samples to calculate H(theta|y,d) samples
		H_theta_posterior = [H_fn(theta) for theta in mcmc_trace]
		
		#Now find the variance of that. Return the inverse as the utility; consider cost separately
		var_H = np.var(H_theta_posterior, ddof=1) #its a sample variance
		U = (1.0/var_H)
		U_trace.append(U)
		print(len(U_trace), end='\r', flush=True)

	Htheta_trace = [H_fn(theta) for theta in mcmc_trace]
	if doPlot:
		plt.xscale('log')
		plt.title("Trace of U")
		plt.plot(U_trace, c='r')
		plt.plot([np.mean(U_trace[0:n]) for n in range(len(U_trace))], c='g')
		plt.plot([np.var(U_trace[0:n], ddof=1) for n in range(len(U_trace))], c='b')
		plt.show()
		
		plt.xscale('log')
		plt.title("Trace of H(theta)")
		plt.plot(Htheta_trace, c='r')
		plt.plot([np.mean(Htheta_trace[0:n]) for n in range(len(Htheta_trace))], c='g')
		plt.plot([np.var(Htheta_trace[0:n], ddof=1) for n in range(len(Htheta_trace))], c='b')
		plt.show()
		
		plt.xscale('log')
		plt.title("Trace of theta")
		plt.plot(mcmc_trace, c='r')
		plt.plot([np.mean(mcmc_trace[0:n]) for n in range(len(mcmc_trace))], c='g')
		plt.plot([np.var(mcmc_trace[0:n], ddof=1) for n in range(len(mcmc_trace))], c='b')
		plt.show()
	
	#Go ahead and save to logfile just in case
	if doLog:
		with open('loop2.csv', 'w+', newline='') as csvfile:
			writer = csv.writer(csvfile)
			write_Utrace = ['']*min + U_trace
			write_Htheta = ['']*min + Htheta_trace
			write_mcmc = mcmc_trace
			for n in range(len(write_mcmc)):
				writer.writerow([n, write_Utrace[n], write_Htheta[n], write_mcmc[n]])


#If you printed from the above function, just run it through here to get plots again
#assumes u = var(H(posterior))
def analyze_trace_loop2(filename):
	U_trace = []
	Htheta_trace = []
	mcmc_trace = []

	print("reading in file...", flush=True)
	with open(filename, 'r', newline='') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if row[1] != '':
				U_trace.append(float(row[1]))
			if row[2] != '':
				Htheta_trace.append(float(row[2]))
			mcmc_trace.append(float(row[3]))

	plt.xscale('log')
	plt.title("Trace of U")
	plt.plot(U_trace, c='r')
	plt.plot([np.mean(U_trace[0:n]) for n in range(len(U_trace))], c='g')
	plt.plot([np.var(U_trace[0:n], ddof=1) for n in range(len(U_trace))], c='b')
	plt.show()
	
	plt.xscale('log')
	plt.title("Trace of H(theta)")
	#plt.plot(Htheta_trace, c='r')
	plt.plot([np.mean(Htheta_trace[0:n]) for n in range(len(Htheta_trace))], c='g')
	plt.plot([np.var(Htheta_trace[0:n], ddof=1) for n in range(len(Htheta_trace))], c='b')
	plt.show()
	
	plt.xscale('log')
	plt.title("Trace of theta")
	#plt.plot(mcmc_trace, c='r')
	plt.plot([np.mean(mcmc_trace[0:n]) for n in range(len(mcmc_trace))], c='g')
	plt.plot([np.var(mcmc_trace[0:n], ddof=1) for n in range(len(mcmc_trace))], c='b')
	plt.show()


#Here, I want to figure out how many likelihood samples are needed to make a good model of likelihood fn
def loop3_convergence(d, exp_fn, H_fn, G_fn, p_theta_rvs, p_theta_pdf, likelihood_est_fn, doLog=True, doPlot=True):
	#likelihood_est_fn should look like calc_likelihood_kernel(d, exp_fn, p_theta_rvs, n1=1000, c='r', showplot=True)
    
	#Come up with a smarter way of determining the fitness of a kde
	0






#####################################################################
# Tools for measuring convergence generally w/ confidence intervals
#####################################################################


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
def is_algorithm_converged(sample, min=1, ci=0.95, closeness=0.95):
	q = scipy.stats.norm.ppf(ci) #phi^-1(ci), quantile fn of x. for example, ci=0.975 means q=1.96
	threshold = 1.0-closeness #this is how many mean-units we want our CI to be
	N = len(sample)
	if N<=min:
		return False
	
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
	0