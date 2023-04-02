#This details an obed script that can be called regardless of the problem
#so it takes as input some function defintions


"""
This function is part of the OBED problem, solving the utility for a given d
This may be an inefficient brute-force approach
This assumes a utility u(d,y,theta) = 1/Var[H(theta|y,d)], minimizing the variance of the HLVA evaluated over the posterior
A higher function must optimize over the outputs of this model to find d*, the d that maximized U_d
d      - the design variable vector to calculate the utility for
exp_fn - function pointer to experiment model, eta(theta, d)
H_fn   - function pointer to high-level verification activity model
G_fn   - function pointer to cost model
"""
def U_brute_varH(d, exp_fn, H_fn, G_fn):
    N1 = 1000 #number of outer loops, ranging over values of y
    
    #Generate a list of y's sampled from likelihood fn, p(y|theta,d)
    #I think you can just do this by running eta(theta,d) over and over - the random epsilon will create the likelihood fn distribution
    #And of course, to get the theta for those executions, you must correspondingly sample from the prior, p(theta)
    
    for y in Y1_list: