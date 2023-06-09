#This details an obed script that can be called regardless of the problem
#so it takes as input some function defintions

#MCMC References:
#https://twiecki.io/blog/2015/11/10/mcmc-sampling/
#https://towardsdatascience.com/bayesian-inference-and-markov-chain-monte-carlo-sampling-in-python-bada1beabca7

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt
import seaborn as sns

"""
This function is part of the OBED problem, solving the utility for a given d
This may be an inefficient brute-force approach
This assumes a utility u(d,y,theta) = 1/Var[H(theta|y,d)], minimizing the variance of the HLVA evaluated over the posterior
A higher function must optimize over the outputs of this model to find d*, the d that maximized U_d
d       - the design variable vector to calculate the utility for
exp_fn  - function pointer to experiment model, eta(theta, d)
H_fn    - function pointer to high-level verification activity model
G_fn    - function pointer to cost model
p_theta_rvs - function pointer to generate samples of the prior distribution p(theta)
p_theta_pdf - function pointer to the pdf of the prior distribution p(theta)
"""
def U_brute_varH(d, exp_fn, H_fn, G_fn, p_theta_rvs, p_theta_pdf, N=1000, burnin=0, lag=1, verbose=False):
    N1 = N #number of outer loops, ranging over values of y
    
    #Generate a list of y's sampled from likelihood fn, p(y|theta,d)
    #I think you can just do this by running eta(theta,d) over and over - the random epsilon will create the likelihood fn distribution
    #And of course, to get the theta for those executions, you must correspondingly sample from the prior, p(theta)
    #so you get a N1-long list of values [y, theta]
    thetas = p_theta_rvs(N1)
    Y1_list = [exp_fn(theta, d) for theta in thetas]
    
    #this outer loop amounts to calculating U over the set of all likely possible data y
    #If you have actual data, i think you can skip or simplify this loop to just do mcmc at that known data?
    
    U_list = []
    for k,y in enumerate(Y1_list):
        if verbose:
            print(k, flush=True)
        #Ok so for each value of y, i need to use MCMC to sample from the posterior
        #I think each of those values of y is my data? In each case, the best i can do for data is one
        #evaluation of the likelihood p(y|theta,d) at a value of y. Looping over the p(y|theta|d)p(theta) values of y
        #makes it valid i think. We'll see
        
        N2 = N*lag + burnin #number of samples of the posterior i want, times lag plus burn-in
        #I will probably have to take into account things like burn-in/convergence and lag? idk
        
        theta_current = p_theta_rvs(1)[0]
        mcmc_trace = []
        for i in range(N2): #MCMC loop
            #Propose a new value of theta from the prior
            theta_proposal = p_theta_rvs(1)[0]
            #Compute likelihood of the "data" for both thetas - this isnt actually easy
            #   I dont have a function of theta that returns p(y|theta,d)
            #   Best I can do is run eta(theta,d) a bunch of times, fit a distribution fn to that, and calculate likelihoods
            #   "probability density estimation"
            N3 = N #number of eta calculations I gotta make for good estimates here
            ylist_current = []
            ylist_proposal = []
            for j in range(N3):
                ylist_current.append(exp_fn(theta_current, d))
                ylist_proposal.append(exp_fn(theta_proposal, d)) #so now we have simulated pdfs of y
            likelihood_current = scipy.stats.gaussian_kde(ylist_current).evaluate(y) #using kde to estimate pdf
            likelihood_proposal = scipy.stats.gaussian_kde(ylist_proposal).evaluate(y) #calc probability at the data y
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
            if i > burnin and i%lag == 0:
                mcmc_trace.append(theta_current)
            
        #Now, use my posterior samples to calculate H(theta|y,d) samples
        H_theta_posterior = [H_fn(theta) for theta in mcmc_trace]
        
        #Now find the variance of that. Return the inverse minus the cost as the utility
        var_H = np.var(H_theta_posterior, ddof=1) #its a sample variance
        U = (1.0/var_H) - G_fn(d)
        U_list.append(U)
        
    return np.average(U_list), U_list
        

def calc_likelihood_kernel(d, exp_fn, p_theta_rvs, n1=1000, c='r', showplot=True):
    thetas = p_theta_rvs(n1)
    Y1_list = [exp_fn(theta, d) for theta in thetas]
    
    y_theta_values = np.vstack([Y1_list, thetas])
    likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values)

    #lets plot this real quick
    if showplot:
        sns.kdeplot(Y1_list, thetas, color=c, shade=True, cmap="Reds", shade_lowest=False)
        plt.show()

    return likelihood_kernel, Y1_list


#This is like U_brute_varH, except we're reusing loop 1 samples in loop 3, like Huan & Marzouk
def U_reloop(d, exp_fn, H_fn, G_fn, p_theta_rvs, p_theta_pdf, utility="varH", n1=10000, n2=1000, burnin=0, lag=1, req=0):   
    #Generate a list of y's sampled from likelihood fn, p(y|theta,d)
    #I think you can just do this by running eta(theta,d) over and over - the random epsilon will create the likelihood fn distribution
    #And of course, to get the theta for those executions, you must correspondingly sample from the prior, p(theta)
    #so you get a N1-long list of values [y, theta]
    
    #Then i want to use my y|d,theta that I generated to create an estimated pdf of y as a fn of theta
    #this is a little more complicated than just a pdf of y like i do originally. Can kde do this?
	
    likelihood_kernel, Y1_list = calc_likelihood_kernel(d, exp_fn, p_theta_rvs, n1, showplot=False)
    
    #this outer loop amounts to calculating U over the set of all likely possible data y
    #If you have actual data, i think you can skip or simplify this loop to just do mcmc at that known data?
    
    U_list = []
    for y in Y1_list:
        #Ok so for each value of y, i need to use MCMC to sample from the posterior
        #I think each of those values of y is my data? In each case, the best i can do for data is one
        #evaluation of the likelihood p(y|theta,d) at a value of y. Looping over the p(y|theta|d)p(theta) values of y
        #makes it valid i think. We'll see
        
        N2 = n2*lag + burnin #number of samples of the posterior i want, times lag plus burn-in
        #I will probably have to take into account things like burn-in/convergence and lag? idk
        
        theta_current = p_theta_rvs(1)[0]
        mcmc_trace = []
        for i in range(N2): #MCMC loop
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
            if i > burnin and i%lag == 0:
                mcmc_trace.append(theta_current)
            
        #Now, use my posterior samples to calculate H(theta|y,d) samples
        H_theta_posterior = [H_fn(theta) for theta in mcmc_trace]
        
		U = 0
		if utility = 'varH':
			#Now find the variance of that posterior list. Return the inverse minus the cost as the utility
			var_H = np.var(H_theta_posterior, ddof=1) #its a sample variance
			U = (1.0/var_H) - G_fn(d)
		elif utility = 'HKL':
			#compute the KL divergence of the sample from posterior to prior
			#Kind of a modified D-optimality
			H_theta_prior = [H_fn(theta) for theta in p_theta_rvs(len(H_theta_posterior))]
			summands = [posterior * np.log(posterior/H_theta_prior[i]) for i,posterior in enumerate(H_theta_posterior))
			U = sum(summands) #D_KL
		elif utility = "Prob95":
			#compute a utility min(P(H > req), 0.95)
			#easily calculated with sample probability, just count em up and find the ratio
			#I think this is functionally the same as optimizing for small variance, just doesnt over-optimize
			num_true = 0
			for h in H_theta_posterior:
				num_true += int(h > req)
			prob = num_true / len(H_theta_posterior)
			U = min(prob, 0.95)
		elif utility = "Prob95_Cost":
			#Using a utility min(P(H > req), 0.95) - Cost*1(P(H > req) >= 0.95)
			#I think this will give us a pareto front that satisfies both meeting the robustness criterion and minimizing cost,
			#but not minimizing the cost until we're meeting the requirement? I think?
			num_true = 0
			for h in H_theta_posterior:
				num_true += int(h > req)
			prob = num_true / len(H_theta_posterior)
			
			isRobust = prob >= 0.95
			penalty = int(isRobust)*G_fn(d) #only apply a penalty cost function when we're in the saturation region; creates convexity along 95% space?
			U = min(prob, 0.95) - penalty #does this need a softening factor? something squared? Idk doesnt look effective
		elif utility = "Cost_95Constraint":
			#Make decreasing cost the utility
			#apply a constraint: if p(H_posterior > req) < 0.95, utility is zero
			#only accepting solutions that satisfy the requirements to the desired certainty
			num_true = 0
			for h in H_theta_posterior:
				num_true += int(h > req)
			prob = num_true / len(H_theta_posterior)
			isReqSatisfied = prob >= 0.95
			
			U = int(isReqSatisfied) * (1.0/G_fn(D))
		else:
			print("something weird")
			exit
        
		
		U_list.append(U)
        
    return np.average(U_list), U_list


	
#Using an approach like Huan & Marzouk 2013, but without polynomial chaos
#This uses D-optimality, the familiar KL divergence of the prior and posterior. Hence a totally different algorithmic approach
def U_reloop_D():