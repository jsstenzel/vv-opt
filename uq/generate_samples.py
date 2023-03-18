#TODO: collapse this down to just consider 1 spectrograph
#TODO: modify this to use the same input as sensitivity_model.py

#import argparse
import os
import sys
import shutil
import csv
import fileinput
sys.path.insert(0, "..")

import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.signal as ss

import numpy as np
import itertools
import math
from SALib.sample import sobol #saltelli, fast_sampler
from copy import deepcopy
import scipy.optimize as optimization


#model is a function that can take x as an input, or x_err  N(0,1) as an input that it applies onto a default x
def sobol_generate_samples(n, var_names, model, doWrite=''):
    N = 2 ** n
    dim = len(var_names)

    #is that what i want?? 
    problem = {
        'names': var_names,
        'num_vars': dim,
        'bounds': [[0, 1]]*dim, #for dists = norm, this is mean=0, sigma=1
        'dists': ['norm']*dim
    }
    
    print("Generating param values...", flush=True)
    param_values = sobol.sample(problem, N, calc_second_order=False)

    ###Setup for func calls
    Y = np.zeros([param_values.shape[0]])

    ###do fn calls
    print("Running model", len(param_values), "times...", flush=True)
    for i, param in enumerate(param_values):
        print('...',i,'...',end='\r',flush=True)
        Y[i] = model(var_names, X_err=param) #this means that we're using the assumed default X, and applying N(0,1) error to it
        
    #write
    if doWrite != '':
        #make a csv file where each line has format [x0, ..., xi, ..., xn, Y]
        writelist = [np.append(val_list, yi) for val_list, yi in zip(param_values,Y)]
        
        with open(doWrite + '.csv', 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for row in writelist:
                writer.writerow(row)
    
    #return [[xn,yn] for xn,yn in zip(param_values, Y)]
    return param_values, Y

def __test_model(var_names, X=[1,2,3], X_err=[0,0,0]):
    x1 = X[0] + X_err[0]
    x2 = X[1] + X_err[1]
    x3 = X[2] + X_err[2]
    return x1 + x2*x3

if __name__ == '__main__':  
    output = sobol_generate_samples(8, ['x1','x2','x3'], __test_model, 'testdata')
    print(output)