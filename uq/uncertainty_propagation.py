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

import numpy as np
import statistics
from scipy.stats import normaltest, norm

#datafile is a csv file where each line has format [x0, ..., xi, ..., xn, Y]
def uncertainty_prop_file(datafile, doPlot=False, doPrint=False):
    y = []
    with open(datafile) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            y.append(float(row[-1]))

    #print(y)
    return uncertainty_prop(y, doPlot, doPrint)

#y_j: list of all model evaluation results, assumed to be single number
def uncertainty_prop(y_j, doPlot=False, doPrint=False):
    mean = statistics.mean(y_j)
    stddev = statistics.stdev(y_j, mean) #note! This is sample stddev, not population stddev. Different n-factor in front

    score, pvalue = normaltest(y_j, nan_policy='raise')
    
    if doPrint:
        print("***Uncertainty propagation results")
        print("Sample mean:",mean)
        print("Sample standard deviation:",stddev)
        print("Gaussian fit p-value:",pvalue,"(p-value > 0.05 means its normal)")
    
    if doPlot:
        n, bins, patches = plt.hist(x=y_j, bins='auto', color='#21aad3', alpha=1.0, rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.ylim(ymax=np.ceil(n.max() / 10) * 10 + n.max()*0.05)
        #plt.xticks(rotation=90)
        plt.xlabel("Y")
        plt.ylabel("Frequency (N=" + str(len(y_j)) + ")")
        
        #add vert line
        plt.axvline(mean, 0, color='b')
        
        plt.show()
    
    return mean, stddev
    
if __name__ == '__main__':  
    print("uncertainty_propagation.py: Use me to analyze and plot the results of model samples!")
    uncertainty_prop_file('testdata.csv', True, True)