#This details an obed script that can be called regardless of the problem
#so it takes as input some function defintions

#MCMC References:
#https://twiecki.io/blog/2015/11/10/mcmc-sampling/
#https://towardsdatascience.com/bayesian-inference-and-markov-chain-monte-carlo-sampling-in-python-bada1beabca7

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt
import seaborn as sns


        

def multivar_likelihood_kernel(d, exp_fn, p_theta_rvs, n1, c='r', showplot=True):
    thetas = p_theta_rvs(n1)
    Y1_list = [exp_fn(theta, d) for theta in thetas]
    
    y_theta_values = np.array([np.concatenate([thetas[i],ys[i]]) for i,_ in enumerate(thetas)])
    likelihood_kernel = scipy.stats.gaussian_kde(y_theta_values.T)

    return likelihood_kernel, Y1_list


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
#This is like U_brute_varH, except we're reusing loop 1 samples in loop 3, like Huan & Marzouk
#Here, we're using straightforward Metropolis-Hastings to solve the multivariate MCMC problem
#This calculates the probability of meeting the requirement in distribution - used as constraint for the cost-optimization
def U_probreq(d, problem, threshold=0.95, n1=10000, n2=1000, burnin=0, lag=1, req=0):   
    #Generate a list of y's sampled from likelihood fn, p(y|theta,d)p(theta)
    
    #Then i want to use my y|d,theta that I generated to create an estimated pdf of y as a fn of theta
    #this is a little more complicated than just a pdf of y like i do originally. Can kde do this?
    likelihood_kernel, Y1_list = multivar_likelihood_kernel(d, problem.eta, problem.prior_rvs, n1, showplot=False)
    
    U_list = []
    for y in Y1_list: #MC loop        
        N2 = n2*lag + burnin #number of samples of the posterior i want, times lag plus burn-in
        #I will probably have to take into account things like burn-in/convergence and lag? idk
        
        theta_current = problem.prior_rvs(1)
        mcmc_trace = []
        for i in range(N2): #MCMC loop
            #Propose a new value of theta from the prior
            theta_proposal = problem.prior_rvs(1)
            #Compute likelihood of the "data" for both thetas - reusing the loop-1 likelihood evaluations
			ytheta_current = np.concatenate([theta_current, y])
			ytheta_proposal = np.concatenate([theta_proposal, y])
            likelihood_current = likelihood_kernel(ytheta_current) #using kde to estimate pdf
            likelihood_proposal = likelihood_kernel(ytheta_proposal) #calc probability at the data y
            #Compute acceptance ratio
            prior_current = np.prod(p_theta_pdf(theta_current))
            prior_proposal = np.prod(p_theta_pdf(theta_proposal))
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
		
		#tweak all this, include req in argument, remove threshold as its dealt with in outer loop
		elif utility = "Prob95": 
			#compute a utility min(P(H > req), 0.95)
			#easily calculated with sample probability, just count em up and find the ratio
			#I think this is functionally the same as optimizing for small variance, just doesnt over-optimize
			num_true = 0
			for h in H_theta_posterior:
				num_true += int(h > req)
			prob = num_true / len(H_theta_posterior)
			U = min(prob, 0.95)

        
		
		U_list.append(U)
        
    return np.average(U_list), U_list

